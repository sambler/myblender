# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

# <pep8 compliant>

import os
import subprocess
import sys
import shutil

# get builder name
if len(sys.argv) < 2:
    sys.stderr.write("Not enough arguments, expecting builder name\n")
    sys.exit(1)

builder = sys.argv[1]

# we run from build/ directory
blender_dir = os.path.join('..', 'blender.git')


def parse_header_file(filename, define):
    import re
    regex = re.compile("^#\s*define\s+%s\s+(.*)" % define)
    with open(filename, "r") as file:
        for l in file:
            match = regex.match(l)
            if match:
                return match.group(1)
    return None

if 'cmake' in builder:
    # cmake

    # Some fine-tuning configuration
    blender_dir = os.path.join('..', blender_dir)
    build_dir = os.path.abspath(os.path.join('..', 'build', builder))
    install_dir = os.path.abspath(os.path.join('..', 'install', builder))
    targets = ['blender']

    chroot_name = None  # If not None command will be delegated to that chroot
    cuda_chroot_name = None  # If not None cuda compilationcommand will be delegated to that chroot
    build_cubins = True  # Whether to build Cycles CUDA kernels
    remove_install_dir = False  # Remove installation folder before building
    bits = 64

    # Config file to be used (relative to blender's sources root)
    cmake_config_file = "build_files/cmake/config/blender_full.cmake"
    cmake_player_config_file = None
    cmake_cuda_config_file = None

    # Set build options.
    cmake_options = []
    cmake_extra_options = ['-DCMAKE_BUILD_TYPE:STRING=Release']
    cuda_cmake_options = []

    if builder.startswith('mac'):
        install_dir = None
        # Set up OSX architecture
        if builder.endswith('x86_64_10_6_cmake'):
            cmake_extra_options.append('-DCMAKE_OSX_ARCHITECTURES:STRING=x86_64')

    elif builder.startswith('win'):
        install_dir = None
        if builder.startswith('win64'):
            cmake_options.append(['-G', '"Visual Studio 12 2013 Win64"'])
        elif builder.startswith('win32'):
            bits = 32
            cmake_options.append(['-G', '"Visual Studio 12 2013"'])

    elif builder.startswith('linux'):
        tokens = builder.split("_")
        glibc = tokens[1]
        if glibc == 'glibc219':
            deb_name = "jessie"
        elif glibc == 'glibc211':
            deb_name = "squeeze"
        remove_install_dir = True
        cmake_config_file = "build_files/buildbot/config/blender_linux.cmake"
        cmake_player_config_file = "build_files/buildbot/config/blender_linux_player.cmake"
        if builder.endswith('x86_64_cmake'):
            chroot_name = 'buildbot_' + deb_name + '_x86_64'
            targets = ['player', 'blender']
        elif builder.endswith('i686_cmake'):
            bits = 32
            chroot_name = 'buildbot_' + deb_name + '_i686'
            cuda_chroot_name = 'buildbot_' + deb_name + '_x86_64'
            targets = ['player', 'blender', 'cuda']

    cmake_options.append("-C" + os.path.join(blender_dir, cmake_config_file))

    # Prepare CMake options needed to configure cuda binaries compilation.
    cuda_cmake_options.append("-DWITH_CYCLES_CUDA_BINARIES=%s" % ('ON' if build_cubins else 'OFF'))
    if build_cubins or 'cuda' in targets:
        if bits == 32:
            cuda_cmake_options.append("-DCUDA_64_BIT_DEVICE_CODE=OFF")
        else:
            cuda_cmake_options.append("-DCUDA_64_BIT_DEVICE_CODE=ON")

    # Only modify common cmake options if cuda doesn't require separate target.
    if 'cuda' not in targets:
        cmake_options += cuda_cmake_options

    if install_dir:
        cmake_options.append("-DCMAKE_INSTALL_PREFIX=%s" % (install_dir))

    cmake_options += cmake_extra_options

    # Prepare chroot command prefix if needed
    if chroot_name:
        chroot_prefix = ['schroot', '-c', chroot_name, '--']
    else:
        chroot_prefix = []
    if cuda_chroot_name:
        cuda_chroot_prefix = ['schroot', '-c', cuda_chroot_name, '--']
    else:
        cuda_chroot_prefix = chroot_prefix[:]

    # Make sure no garbage remained from the previous run
    # (only do it if builder requested this)
    if remove_install_dir:
        if os.path.isdir(install_dir):
            shutil.rmtree(install_dir)

    for target in targets:
        print("Building target %s" % (target))
        # Construct build directory name based on the target
        target_build_dir = build_dir
        target_chroot_prefix = chroot_prefix[:]
        if target != 'blender':
            target_build_dir += '_' + target
        target_name = 'install'
        # Make sure build directory exists and enter it
        if not os.path.isdir(target_build_dir):
            os.mkdir(target_build_dir)
        os.chdir(target_build_dir)
        # Tweaking CMake options to respect the target
        target_cmake_options = cmake_options[:]
        if target == 'player':
            target_cmake_options.append("-C" + os.path.join(blender_dir, cmake_player_config_file))
        elif target == 'cuda':
            target_cmake_options += cuda_cmake_options
            target_chroot_prefix = cuda_chroot_prefix[:]
            target_name = 'cycles_kernel_cuda'
        # If cuda binaries are compiled as a separate target, make sure
        # other targets don't compile cuda binaries.
        if 'cuda' in targets and target != 'cuda':
            target_cmake_options.append("-DWITH_CYCLES_CUDA_BINARIES=OFF")
        # Configure the build
        print("CMake options:")
        print(target_cmake_options)
        if os.path.exists('CMakeCache.txt'):
            print("Removing CMake cache")
            os.remove('CMakeCache.txt')
        retcode = subprocess.call(target_chroot_prefix + ['cmake', blender_dir] + target_cmake_options)
        if retcode != 0:
            print('Condifuration FAILED!')
            sys.exit(retcode)

        if 'win32' in builder:
            command = ['msbuild', 'INSTALL.vcxproj', '/Property:PlatformToolset=v120_xp', '/p:Configuration=Release']
        elif 'win64' in builder:
            command = ['msbuild', 'INSTALL.vcxproj', '/p:Configuration=Release']
        else:
            command = target_chroot_prefix + ['make', '-s', '-j2', target_name]

        print("Executing command:")
        print(command)
        retcode = subprocess.call(command)

        if retcode != 0:
            sys.exit(retcode)

        if builder.startswith('linux') and target == 'cuda':
            blender_h = os.path.join(blender_dir, "source", "blender", "blenkernel", "BKE_blender.h")
            blender_version = int(parse_header_file(blender_h, 'BLENDER_VERSION'))
            blender_version = "%d.%d" % (blender_version // 100, blender_version % 100)
            kernels = os.path.join(target_build_dir, 'intern', 'cycles', 'kernel')
            install_kernels = os.path.join(install_dir, blender_version, 'scripts', 'addons', 'cycles', 'lib')
            os.mkdir(install_kernels)
            print("Copying cuda binaries from %s to %s" % (kernels, install_kernels))
            os.system('cp %s/*.cubin %s' % (kernels, install_kernels))

else:
    print("Unknown building system")
    sys.exit(1)
