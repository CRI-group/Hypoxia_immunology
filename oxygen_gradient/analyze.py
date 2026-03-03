import os
import time
import pickle
import functools
import subprocess
import xml.etree.ElementTree as ET
import pandas as pd
import cv2
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from tqdm import tqdm
from pprint import pprint

from argparse import ArgumentParser, ArgumentTypeError


def timeit(func):
    """Decorator to time function execution"""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} executed in {end_time - start_time:.4f} seconds")
        return result

    return wrapper


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise ArgumentTypeError("Boolean value expected.")


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument("path_xml", help="Path to xml", type=str)
    parser.add_argument("path_imgs_dir", help="Path to images dir", type=str)
    args = vars(parser.parse_args())
    return args


def check_arguments(args):
    pass


class WSITrackViewer:
    def __init__(self, xml_path, image_dir, fade_range=4, cache_file="xml_cache.pkl"):
        self.xml_path = xml_path
        self.image_dir = image_dir
        self.fade_range = fade_range
        self.cache_file = cache_file

        # Load image files sorted by time point
        self.image_files = []
        for f in os.listdir(image_dir):
            if f.endswith((".tif", ".png", ".jpg")):
                self.image_files.append(os.path.join(image_dir, f))
        self.image_files.sort()
        self.num_frames = len(self.image_files)

        # get scale
        self.scale = self.parse_scale(xml_path)
        # get image size
        self.img_height, self.img_width = self.load_image(0).shape

        # Load and parse XML
        df_spots, df_tracks, df_edges = self.parse_xml()
        df_edges = self.insert_x_and_y_coords_into_df_edges(df_edges, df_spots)
        self.df_spots = df_spots
        self.df_tracks = df_tracks
        self.df_edges = df_edges

        # # visualize spots on a given frame
        # self.visualize_spots_at_frame(0)

        # ANALYSIS

        # 1 measure mean shift
        self.measure_mean_shift(df_edges)

        # # 2 visualize track ends
        # track_ends_output_dir = f"{image_dir.rstrip('/')}_track_ends"
        # self.visualize_track_ends(df_edges, df_spots, self.image_files, track_ends_output_dir)

        # 3 count cell counts at diffeent segments for the final frame/measurement
        path_out_vis = os.path.join(f"{image_dir.rstrip('/')}_segment_visualization.png")
        df_spots = self.measure_segment_counts_at_final_frame(
            df_spots,
            self.img_width,
            self.img_height,
            N_segments=5,
            left_margin=0.08,
            right_margin=0.06,
            rotation=6.0,
            path_visualization_image_out=path_out_vis,
        )

        # # Set up visualization
        # self.current_frame = 0
        # self.fig, self.ax = plt.subplots()
        # plt.subplots_adjust(bottom=0.2)
        # # Load first image
        # self.img = self.load_image(self.current_frame)
        # self.img_display = self.ax.imshow(self.img, cmap="gray")
        # self.overlay_tracks()
        # # Add slider
        # self.slider_ax = plt.axes([0.2, 0.05, 0.65, 0.03])
        # self.slider = Slider(self.slider_ax, "Frame", 0, self.num_frames - 1, valinit=0, valfmt="%0.0f")
        # self.slider.on_changed(self.update)
        # plt.show()

    def measure_mean_shift(self, df_edges):
        df_edges["X_diff"] = df_edges.X_target - df_edges.X_source
        df_edges["Y_diff"] = df_edges.Y_target - df_edges.Y_source

        # compute average shift per track
        df = df_edges.groupby("Track_ID").agg("mean")
        # compute average shift over the whole stack
        average_shift_in_px = df[["X_diff", "Y_diff"]].mean() / self.scale

        print("")
        print("Average shift of tracked cells in pixels over the image stack")
        print(pd.DataFrame(average_shift_in_px).T)

    def visualize_track_ends(self, df_edges, df_spots, image_files, output_dir):
        os.makedirs(output_dir, exist_ok=True)

        last_frame = df_spots.Frame.max()
        # get edge ends
        df = df_edges.groupby("Track_ID").max()[["X_target", "Y_target", "Frame_target"]]
        # get track ends that are not on the last frame
        track_ends = df[df.Frame_target != last_frame]

        # def close_on_click(event):
        #     plt.close()

        for frame_idx, fn in enumerate(tqdm(image_files, "saving track end images")):

            # fig, ax = plt.subplots()
            # fig.canvas.mpl_connect("key_press_event", close_on_click)
            fn_base = os.path.basename(fn)
            fn_base = ".".join(os.path.splitext(fn_base)[:-1])

            img = self.load_image(frame_idx)

            fig, ax = plt.subplots()
            ax.imshow(img, cmap="gray")
            df_frame = track_ends[track_ends.Frame_target == frame_idx]
            ax.scatter(
                df_frame.X_target / self.scale,
                df_frame.Y_target / self.scale,
                color="red",
                marker="x",
                s=7,
                linewidth=0.5,
            )
            # Remove all axes and borders
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_frame_on(False)
            path_out = os.path.join(output_dir, f"frame_{frame_idx}_{fn_base}.png")
            plt.savefig(path_out, dpi=300, bbox_inches="tight", pad_inches=0)
            plt.close(fig)

    def measure_segment_counts_at_final_frame(
        self,
        df_spots,
        img_width,
        img_height,
        N_segments,
        left_margin,
        right_margin,
        rotation,
        path_visualization_image_out,
    ):
        """
        Divide image vertically into segments, starting from separate left and right margins,
        include rotation, and measure cell counts at the final frame.
        """

        # Define separate margins
        left_margin = img_width * left_margin
        right_margin = img_width * (1 - right_margin)

        # Compute segment boundaries **between** the margins
        segment_width = (right_margin - left_margin) / N_segments
        segment_edges = [left_margin + i * segment_width for i in range(N_segments + 1)]

        def rotate_point(x, y, degrees):
            """Apply rotation transformation to a point (x, y)"""
            theta = np.radians(degrees)
            cos_t, sin_t = np.cos(theta), np.sin(theta)

            x_center, y_center = img_width / 2, img_height / 2  # Rotate around image center
            x_shifted, y_shifted = x - x_center, y - y_center
            x_rot = x_shifted * cos_t - y_shifted * sin_t + x_center
            y_rot = x_shifted * sin_t + y_shifted * cos_t + y_center
            return x_rot, y_rot

        # Rotate segment boundaries to match rotated space
        rotated_segment_edges = [
            rotate_point(x, img_height / 2, rotation)[0] for x in segment_edges
        ]

        # Visualize segments on the last image
        img = self.load_image(df_spots.Frame.max())
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(img, cmap="gray")

        # Draw vertical segment boundaries in red (with rotation)
        for i in range(1, N_segments):
            x_pos = segment_edges[i]
            x1, y1 = rotate_point(x_pos, 0, rotation)  # Top of the image
            x2, y2 = rotate_point(x_pos, img_height, rotation)  # Bottom of the image
            ax.plot([x1, x2], [y1, y2], color="red", linestyle="--", linewidth=1, label=None if i > 1 else "Segments")

        # Draw left and right margins in green (with rotation)
        for margin, label in zip([left_margin, right_margin], ["Margins", None]):
            x1, y1 = rotate_point(margin, 0, rotation)
            x2, y2 = rotate_point(margin, img_height, rotation)
            ax.plot([x1, x2], [y1, y2], color="green", linestyle="--", linewidth=1, label=label)

        ax.legend(loc="upper right")
        ax.set_title(f"Segment Visualization (Rotation={rotation}°)")
        ax.set_axis_off()

        # Function to assign segments considering rotation
        def assign_segment(x, y):
            """Assign points to segments based on rotated boundaries"""
            x_rot, _ = rotate_point(x, y, rotation)  # Rotate point FORWARD to match segment boundaries

            if x_rot < rotated_segment_edges[0]:  # Compare in rotated space
                return 1  # Points left of left margin belong to first segment
            if x_rot >= rotated_segment_edges[-1]:
                return N_segments  # Points right of right margin belong to last segment

            for i in range(len(rotated_segment_edges) - 1):
                if rotated_segment_edges[i] <= x_rot < rotated_segment_edges[i + 1]:
                    return i + 1  # Segments numbered from 1 to N_segments

            return None  # In case X is out of range

        # Assigning segments
        df_spots["Segment"] = df_spots.apply(lambda row: assign_segment(row["X"], row["Y"]), axis=1)

        last_frame = df_spots.Frame.max()
        last_frame_spots = df_spots[df_spots.Frame == last_frame]
        spots_in_segment = last_frame_spots[last_frame_spots.Segment == 1]
        ax.scatter(
            spots_in_segment.X / self.scale,
            spots_in_segment.Y / self.scale,
            color="red",
            marker="x",
            s=7,
            linewidth=0.5,
        )

        plt.savefig(path_visualization_image_out, dpi=300)
        plt.close()

    def measure_segment_counts_at_final_frame(
        self,
        df_spots,
        img_width,
        img_height,
        N_segments,
        left_margin,
        right_margin,
        rotation,
        path_visualization_image_out,
    ):
        """
        Divide image vertically into segments, starting from separate left and right margins,
        include rotation, and measure cell counts at the final frame.
        """
        
        # Define separate margins
        # image size is already in the correct scale, no need to use self.scale
        # only points are in the subsampled scale
        left_margin = img_width * left_margin
        right_margin = img_width * (1 - right_margin)

        # Compute segment boundaries **between** the margins
        segment_width = (right_margin - left_margin) / N_segments
        segment_edges = [left_margin + i * segment_width for i in range(N_segments + 1)]

        def rotate_point(x, y, degrees):
            """Apply rotation transformation to a point (x, y)"""

            # Rotation matrix
            theta = np.radians(degrees)
            cos_t, sin_t = np.cos(theta), np.sin(theta)

            x_center, y_center = img_width / 2, img_height / 2  # Rotate around image center
            x_shifted, y_shifted = x - x_center, y - y_center
            x_rot = x_shifted * cos_t - y_shifted * sin_t + x_center
            y_rot = x_shifted * sin_t + y_shifted * cos_t + y_center
            return x_rot, y_rot

        # Visualize segments on the last image
        img = self.load_image(df_spots.Frame.max())
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(img, cmap="gray")
        # Draw vertical segment boundaries in red (with rotation)
        for i in range(1, N_segments):
            x_pos = segment_edges[i]
            x1, y1 = rotate_point(x_pos, 0, rotation)  # Top of the image
            x2, y2 = rotate_point(x_pos, img_height, rotation)  # Bottom of the image
            ax.plot([x1, x2], [y1, y2], color="red", linestyle="--", linewidth=1, label=None if i > 1 else f"Segments")
        # Draw left and right margins in green (with rotation)
        for margin, label in zip([left_margin, right_margin], ["Margins", None]):
            x1, y1 = rotate_point(margin, 0, rotation)
            x2, y2 = rotate_point(margin, img_height, rotation)
            ax.plot([x1, x2], [y1, y2], color="green", linestyle="--", linewidth=1, label=label)
        ax.legend(loc="upper right")
        ax.set_title(f"Segment Visualization (Rotation={rotation}°)")
        ax.set_axis_off()
        # plt.show()

        # Function to assign segments considering rotation
        def assign_segment(x, y):
            """Assign points to segments based on rotated boundaries"""
            x_rot, _ = rotate_point(x / self.scale, y / self.scale, -rotation)  # Rotate point to match rotated segment positions

            if x_rot < left_margin:
                return 1  # Points left of left margin belong to first segment
            if x_rot >= right_margin:
                return N_segments  # Points right of right margin belong to last segment

            for i in range(len(segment_edges) - 1):
                if segment_edges[i] <= x_rot < segment_edges[i + 1]:
                    return i + 1  # Segments numbered from 1 to N_segments

            return None  # In case X is out of range

        # Assigning segments
        df_spots["Segment"] = df_spots.apply(lambda row: assign_segment(row["X"], row["Y"]), axis=1)

        # # visualize points on a specific segment
        # last_frame = df_spots.Frame.max()
        # last_frame_spots = df_spots[df_spots.Frame == last_frame]
        # spots_in_segment = last_frame_spots[last_frame_spots.Segment == 1]
        # ax.scatter(spots_in_segment.X / self.scale, spots_in_segment.Y / self.scale, color="red", marker="x", s=7, linewidth=0.5)

        plt.savefig(path_visualization_image_out, dpi=300)
        plt.close()

        cells_at_first = pd.DataFrame(df_spots[df_spots.Frame == df_spots.Frame.min()].groupby("Segment").count()["ID"]).rename(columns={"ID":"cell count"}).T
        cells_at_last = pd.DataFrame(df_spots[df_spots.Frame == df_spots.Frame.max()].groupby("Segment").count()["ID"]).rename(columns={"ID":"cell count"}).T

        print("Cell count in each segment at the first measurement")
        print(cells_at_first)

        print("Cell count in each segment at the last measurement")
        print(cells_at_last)

        return df_spots


    def visualize_spots_at_frame(self, frame_idx):
        img = self.load_image(frame_idx)
        plt.imshow(img)
        f0_spots = self.df_spots[self.df_spots.Frame == 0]
        plt.scatter(f0_spots.X / self.scale, f0_spots.Y / self.scale)
        plt.show()

    def insert_x_and_y_coords_into_df_edges(self, df_edges, df_spots):
        # Merge for source spots
        df_edges = df_edges.merge(
            df_spots[["ID", "X", "Y", "Frame"]], left_on="SPOT_SOURCE_ID", right_on="ID", how="left"
        )
        df_edges.rename(columns={"X": "X_source", "Y": "Y_source", "Frame": "Frame_source"}, inplace=True)
        df_edges.drop(columns=["ID"], inplace=True)  # Drop redundant 'ID' column

        # Merge for target spots
        df_edges = df_edges.merge(
            df_spots[["ID", "X", "Y", "Frame"]], left_on="SPOT_TARGET_ID", right_on="ID", how="left"
        )
        df_edges.rename(columns={"X": "X_target", "Y": "Y_target", "Frame": "Frame_target"}, inplace=True)
        df_edges.drop(columns=["ID"], inplace=True)  # Drop redundant 'ID' column

        return df_edges

    def parse_scale(self, xml_path):
        # Parse scale from xml
        print("Parse data scale from xml file")
        cmd = f'grep -Po "(?<=dx = )([0-9\.]*)" {xml_path}'
        dx = subprocess.check_output(cmd, shell=True, text=True)
        cmd = f'grep -Po "(?<=dy = )([0-9\.]*)" {xml_path}'
        dy = subprocess.check_output(cmd, shell=True, text=True)
        cmd = f'grep -Po "(?<=dy = )([0-9\.]*)" {xml_path}'
        dz = subprocess.check_output(cmd, shell=True, text=True)
        assert dx == dy == dz, "Scale mismatch for different image dimensions"
        scale = float(dx)
        return scale

    @timeit
    def load_or_parse_xml(self):
        """Loads the cached XML parsing result or parses the XML file if cache is unavailable."""
        if os.path.exists(self.cache_file):
            print("loading from cache")
            with open(self.cache_file, "rb") as f:
                return pickle.load(f)
        else:
            print("No cache found, parsing from xml..")
            df = self.parse_xml()
            with open(self.cache_file, "wb") as f:
                pickle.dump(df, f)
            print("parsing done")
            return df

    @timeit
    def parse_xml(self):
        """Parses the TrackMate XML file into a DataFrame."""
        tree = ET.parse(self.xml_path)
        root = tree.getroot()

        spots = []
        for frame in root.findall(".//SpotsInFrame"):
            frame_number = int(frame.get("frame"))
            for spot in frame.findall("Spot"):
                spot_data = {
                    "ID": int(spot.get("ID")),
                    "name": spot.get("ID"),
                    "Frame": frame_number,
                    "QUALITY": float(spot.get("QUALITY")),
                    "X": float(spot.get("POSITION_X")),
                    "Y": float(spot.get("POSITION_Y")),
                    "Z": float(spot.get("POSITION_Z")),
                    "T": float(spot.get("POSITION_T")),
                    "MIN_INTENSITY_CH1": float(spot.get("MIN_INTENSITY_CH1")),
                    "TOTAL_INTENSITY_CH1": float(spot.get("TOTAL_INTENSITY_CH1")),
                    "CONTRAST_CH1": float(spot.get("CONTRAST_CH1")),
                    "SNR_CH1": float(spot.get("SNR_CH1")),
                    "FRAME": int(spot.get("FRAME")),
                    "MEDIAN_INTENSITY_CH1": float(spot.get("MEDIAN_INTENSITY_CH1")),
                    "VISIBILITY": int(spot.get("VISIBILITY")),
                    "RADIUS": float(spot.get("RADIUS")),
                    "MEAN_INTENSITY_CH1": float(spot.get("MEAN_INTENSITY_CH1")),
                    "MAX_INTENSITY_CH1": float(spot.get("MAX_INTENSITY_CH1")),
                    "STD_INTENSITY_CH1": float(spot.get("STD_INTENSITY_CH1")),
                }
                spots.append(spot_data)
        df_spots = pd.DataFrame(spots)

        tracks = []
        for track in root.findall(".//Track"):
            track_id = int(track.get("TRACK_ID"))
            track_data = {
                "TRACK_ID": int(track.get("TRACK_ID")),
                "TRACK_INDEX": int(track.get("TRACK_INDEX")),
                "NUMBER_SPOTS": int(track.get("NUMBER_SPOTS")),
                "NUMBER_GAPS": int(track.get("NUMBER_GAPS")),
                "NUMBER_SPLITS": int(track.get("NUMBER_SPLITS")),
                "NUMBER_MERGES": int(track.get("NUMBER_MERGES")),
                "NUMBER_COMPLEX": int(track.get("NUMBER_COMPLEX")),
                "LONGEST_GAP": int(track.get("LONGEST_GAP")),
                "TRACK_DURATION": float(track.get("TRACK_DURATION")),
                "TRACK_START": float(track.get("TRACK_START")),
                "TRACK_STOP": float(track.get("TRACK_STOP")),
                "TRACK_DISPLACEMENT": float(track.get("TRACK_DISPLACEMENT")),
                "TRACK_X_LOCATION": float(track.get("TRACK_X_LOCATION")),
                "TRACK_Y_LOCATION": float(track.get("TRACK_Y_LOCATION")),
                "TRACK_Z_LOCATION": float(track.get("TRACK_Z_LOCATION")),
                "TRACK_MEAN_SPEED": float(track.get("TRACK_MEAN_SPEED")),
                "TRACK_MAX_SPEED": float(track.get("TRACK_MAX_SPEED")),
                "TRACK_MIN_SPEED": float(track.get("TRACK_MIN_SPEED")),
                "TRACK_MEDIAN_SPEED": float(track.get("TRACK_MEDIAN_SPEED")),
                "TRACK_STD_SPEED": float(track.get("TRACK_STD_SPEED")),
                "TRACK_MEAN_QUALITY": float(track.get("TRACK_MEAN_QUALITY")),
                "TOTAL_DISTANCE_TRAVELED": float(track.get("TOTAL_DISTANCE_TRAVELED")),
                "MAX_DISTANCE_TRAVELED": float(track.get("MAX_DISTANCE_TRAVELED")),
                "CONFINEMENT_RATIO": float(track.get("CONFINEMENT_RATIO")),
                "MEAN_STRAIGHT_LINE_SPEED": float(track.get("MEAN_STRAIGHT_LINE_SPEED")),
                "LINEARITY_OF_FORWARD_PROGRESSION": float(track.get("LINEARITY_OF_FORWARD_PROGRESSION")),
                "MEAN_DIRECTIONAL_CHANGE_RATE": float(track.get("MEAN_DIRECTIONAL_CHANGE_RATE")),
            }
            tracks.append(track_data)
        df_tracks = pd.DataFrame(tracks)

        edges = []
        for track in root.findall(".//Track"):
            track_id = int(track.get("TRACK_ID"))
            for edge in track.findall("Edge"):
                edge_data = {
                    "Track_ID": track_id,
                    "SPOT_SOURCE_ID": int(edge.get("SPOT_SOURCE_ID")),
                    "SPOT_TARGET_ID": int(edge.get("SPOT_TARGET_ID")),
                    "SPEED": float(edge.get("SPEED")),
                    "LINK_COST": float(edge.get("LINK_COST")),
                    "DIRECTIONAL_CHANGE_RATE": float(edge.get("DIRECTIONAL_CHANGE_RATE")),
                    "DISPLACEMENT": float(edge.get("DISPLACEMENT")),
                    "EDGE_TIME": float(edge.get("EDGE_TIME")),
                    "EDGE_X_LOCATION": float(edge.get("EDGE_X_LOCATION")),
                    "EDGE_Y_LOCATION": float(edge.get("EDGE_Y_LOCATION")),
                    "EDGE_Z_LOCATION": float(edge.get("EDGE_Z_LOCATION")),
                }
                edges.append(edge_data)
        df_edges = pd.DataFrame(edges)

        return df_spots, df_tracks, df_edges

    def load_image(self, frame_idx):
        """Loads an image efficiently."""
        img = cv2.imread(self.image_files[frame_idx], cv2.IMREAD_GRAYSCALE)
        return img

    def overlay_tracks(self):
        """Overlays tracks on the current image."""
        self.ax.clear()
        self.img_display = self.ax.imshow(self.img, cmap="gray")

        for offset in tqdm(range(-self.fade_range, self.fade_range + 1), "drawing tracks"):
            faded_frame = self.current_frame + offset
            if 0 <= faded_frame < self.num_frames:
                alpha = 1.0 - abs(offset) / (self.fade_range + 1)
                frame_tracks = self.df[self.df["Frame"] == faded_frame]
                for _, row in frame_tracks.iterrows():
                    self.ax.scatter(row["X"], row["Y"], c="red", alpha=alpha, s=10)

        self.ax.set_title(f"Frame {self.current_frame}")
        self.fig.canvas.draw()

    def update(self, val):
        """Updates the visualization when slider is changed."""
        self.current_frame = int(self.slider.val)
        self.img = self.load_image(self.current_frame)
        self.overlay_tracks()


def main():
    args = parse_arguments()

    xml_path = args["path_xml"]
    path_cache = f"{os.path.basename(xml_path)}_cache.pkl"
    image_dir = args["path_imgs_dir"]
    viewer = WSITrackViewer(xml_path, image_dir, cache_file=path_cache)


if __name__ == "__main__":
    main()
