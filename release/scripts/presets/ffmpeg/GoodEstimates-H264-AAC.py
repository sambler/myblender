import bpy

bpy.context.scene.render.ffmpeg.format = "H264"
#bpy.context.scene.render.ffmpeg.codec = "H264"
bpy.context.scene.render.ffmpeg.use_lossless_output = True

bpy.context.scene.render.ffmpeg.gopsize = bpy.context.scene.render.fps / 2

if resolution.x == 4096:
    bpy.context.scene.render.ffmpeg.video_bitrate = 40000 * 1.2
    bpy.context.scene.render.ffmpeg.maxrate = 40000 * 1.5
    bpy.context.scene.render.ffmpeg.minrate = 0
else if resolution.x == 2560:
    bpy.context.scene.render.ffmpeg.video_bitrate = 10000 * 1.2
    bpy.context.scene.render.ffmpeg.maxrate = 10000 * 1.5
    bpy.context.scene.render.ffmpeg.minrate = 0
else if resolution.x == 1920:
    bpy.context.scene.render.ffmpeg.video_bitrate = 8000 * 1.2
    bpy.context.scene.render.ffmpeg.maxrate = 8000 * 1.5
    bpy.context.scene.render.ffmpeg.minrate = 0
else if resolution.x == 1280:
    bpy.context.scene.render.ffmpeg.video_bitrate = 5000 * 1.2
    bpy.context.scene.render.ffmpeg.maxrate = 5000 * 1.5
    bpy.context.scene.render.ffmpeg.minrate = 0
else if resolution.x == 640:
    bpy.context.scene.render.ffmpeg.video_bitrate = 2500 * 1.2
    bpy.context.scene.render.ffmpeg.maxrate = 2500 * 1.5
    bpy.context.scene.render.ffmpeg.minrate = 0
else if resolution.x == 480:
    bpy.context.scene.render.ffmpeg.video_bitrate = 1000 * 1.2
    bpy.context.scene.render.ffmpeg.maxrate = 1000 * 1.5
    bpy.context.scene.render.ffmpeg.minrate = 0


bpy.context.scene.render.ffmpeg.buffersize = 224 * 8
bpy.context.scene.render.ffmpeg.packetsize = 2048
bpy.context.scene.render.ffmpeg.muxrate = 10080000

bpy.context.scene.render.ffmpeg.audio_codec = "AAC"
bpy.context.scene.render.ffmpeg.audio_bitrate = 384
bpy.context.scene.render.ffmpeg.audio_mixrate = 48000
bpy.context.scene.render.ffmpeg.audio_channels = "STEREO"
