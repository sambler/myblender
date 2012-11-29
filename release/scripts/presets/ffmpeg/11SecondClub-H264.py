import bpy

bpy.context.scene.render.ffmpeg.format = "H264"
#bpy.context.scene.render.ffmpeg.codec = "H264"
bpy.context.scene.render.ffmpeg.use_lossless_output = True

bpy.context.scene.render.ffmpeg.gopsize = 18

bpy.context.scene.render.ffmpeg.video_bitrate = 6000
bpy.context.scene.render.ffmpeg.maxrate = 9000
bpy.context.scene.render.ffmpeg.minrate = 0
bpy.context.scene.render.ffmpeg.buffersize = 224 * 8
bpy.context.scene.render.ffmpeg.packetsize = 2048
bpy.context.scene.render.ffmpeg.muxrate = 10080000

bpy.context.scene.render.ffmpeg.audio_codec = "MP3"
bpy.context.scene.render.ffmpeg.audio_bitrate = 128
bpy.context.scene.render.ffmpeg.audio_mixrate = 48000
bpy.context.scene.render.ffmpeg.audio_channels = "STEREO"

