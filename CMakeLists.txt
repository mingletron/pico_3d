cmake_minimum_required(VERSION 3.13)
include(pico_sdk_import.cmake)
project(draw3d VERSION 1.0.0)
add_executable(draw3d
               draw3d.cpp
               ssd1306.cpp)

pico_sdk_init()

pico_enable_stdio_usb(draw3d 1)
pico_enable_stdio_uart(draw3d 1)
pico_add_extra_outputs(draw3d)

target_link_libraries(draw3d
                      pico_stdlib
                      hardware_gpio
                      hardware_i2c
                      hardware_adc
                      pico_binary_info)
