from valis import registration, feature_detectors, feature_matcher
from valis.micro_rigid_registrar import MicroRigidRegistrar

micro_rigid_registrar_cls = MicroRigidRegistrar

micro_rigid_registrar_params = {
    "scale": 0.12,
    "tile_wh": 1024,
    "roi": "mask"
}

slide_src_dir = "/mnt/project/DAPIs"
results_dst_dir = "/mnt/project/VALIS/valis_files/micro/"
registered_slide_dst_dir = "/mnt/project/VALIS/registered_DAPIs_micro"
reference_slide = "/mnt/project/DAPIs/Astro_15B_CA9_CD45_DAPI_x20_z0.tif"

registrar = registration.Valis(
    slide_src_dir,
    results_dst_dir,
    reference_img_f=reference_slide,
    max_processed_image_dim_px=2048,
    align_to_reference=True,
    micro_rigid_registrar_cls=micro_rigid_registrar_cls,
    micro_rigid_registrar_params=micro_rigid_registrar_params
)

rigid_registrar, non_rigid_registrar, error_df = registrar.register()

# not needed. This does non-rigid micro registration.
# rigid micro done automatically inside registrar.register()
# registrar.register_micro()

registrar.warp_and_save_slides(
    registered_slide_dst_dir,
    non_rigid=False
)

registration.kill_jvm()
