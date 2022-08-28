mkdir -p bin/debug/
cp -u shaders bin/debug/ -d -r
cp -u textures bin/debug/ -d -r
gcc -o bin/debug/hazelholm_c stb_image.c delo2d.c mpl.c mpl_debug.c game.c main.c -Ilibs -lglfw -lGLEW -lGL -msse2 -lm -ldl -O2
./bin/debug/hazelholm_c