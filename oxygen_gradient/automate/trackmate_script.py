import sys
import os
import java.lang.Runtime  # for increasing memory at runtime
import java.lang.System  # for increasing memory at runtime
from ij import IJ
from ij.plugin import FolderOpener
from fiji.plugin.trackmate import Model, Settings, TrackMate, Logger
from fiji.plugin.trackmate.io import TmXmlWriter
from fiji.plugin.trackmate.detection import LogDetectorFactory
from fiji.plugin.trackmate.tracking.jaqaman import SparseLAPTrackerFactory
# from fiji.plugin.trackmate.features.FeatureFilter import FeatureFilter
from ij import WindowManager

# Increase memory allocation before image processing
java.lang.System.setProperty('java.vm.heap.max', '14G')
java.lang.System.setProperty('java.vm.heap.min', '2G')
# Print to check max heap size
print("Max Heap Size (MB):", java.lang.Runtime.getRuntime().maxMemory() / (1024*1024))


# Reload UTF-8 settings to avoid encoding issues
reload(sys)
sys.setdefaultencoding('utf-8')

# Open image sequence from the specified folder
image_path = "/home/masi/Projects/oxygen_gradient/data/test_sample/microglia_raw_images-acidic_region1/"
options = "sort use virtual"
imp = FolderOpener.open(image_path, options)

if imp is None:
    print("Error: Could not open image sequence from path:", image_path)
    sys.exit(1)

# Initialize model
model = Model()
model.setLogger(Logger.IJ_LOGGER)



# Configure detector settings
settings = Settings(imp)
settings.detectorFactory = LogDetectorFactory()
settings.detectorSettings = {
    'DO_SUBPIXEL_LOCALIZATION': True,
    'DO_MEDIAN_FILTERING': False,
    'RADIUS': 17.0,
    'THRESHOLD': 73.0,
    'TARGET_CHANNEL': 1,
}


# # Apply spot filter
# filter1 = FeatureFilter('QUALITY', 300, True)
# settings.addSpotFilter(filter1)

# Configure tracking settings
settings.trackerFactory = SparseLAPTrackerFactory()
settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
settings.trackerSettings['ALLOW_TRACK_SPLITTING'] = True
settings.trackerSettings['ALLOW_TRACK_MERGING'] = True
settings.trackerSettings['ALLOW_GAP_CLOSING'] = False

settings.addAllAnalyzers()

# Instantiate TrackMate plugin
trackmate = TrackMate(model, settings)

# Run tracking
if not trackmate.checkInput():
    sys.exit(trackmate.getErrorMessage())
if not trackmate.process():
    sys.exit(trackmate.getErrorMessage())

# Define output XML path
output_xml_path = os.path.join(os.getcwd(), "trackmate_output.xml")

# Write results to XML
writer = TmXmlWriter(output_xml_path)
writer.appendModel(model)
writer.writeToFile()

print("Saved track data to XML: " + output_xml_path)
