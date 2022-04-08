cmake_minimum_required (VERSION 3.14)

project(CMSISBUILD)

# DSP Sources
SET(DSP ${ROOT}/DSP)

list(APPEND CMAKE_MODULE_PATH ${DSP}/Include)
list(APPEND CMAKE_MODULE_PATH ${DSP})

include(arm_math.h)

