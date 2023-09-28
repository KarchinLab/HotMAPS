set(PULSEAUDIO_FOUND TRUE)

set(PULSEAUDIO_VERSION_MAJOR 16)
set(PULSEAUDIO_VERSION_MINOR 1)
set(PULSEAUDIO_VERSION 16.1)
set(PULSEAUDIO_VERSION_STRING "16.1")

find_path(PULSEAUDIO_INCLUDE_DIR pulse/pulseaudio.h HINTS "/lab/projects1/GRUMP/HotMAPS2/.snakemake/conda/8753e9016cc1f5b96951553cc62b9320_/include")
find_library(PULSEAUDIO_LIBRARY NAMES pulse libpulse HINTS "/lab/projects1/GRUMP/HotMAPS2/.snakemake/conda/8753e9016cc1f5b96951553cc62b9320_/lib")
find_library(PULSEAUDIO_MAINLOOP_LIBRARY NAMES pulse-mainloop-glib libpulse-mainloop-glib HINTS "/lab/projects1/GRUMP/HotMAPS2/.snakemake/conda/8753e9016cc1f5b96951553cc62b9320_/lib")
