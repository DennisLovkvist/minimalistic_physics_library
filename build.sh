mkdir -p bin/debug/
cp -u shaders bin/debug/ -d -r
cp -u textures bin/debug/ -d -r
gcc -o bin/debug/hazelholm_c stb_image.c delo2d.c delo_math.c mpl.c main.c -Ilibs -lglfw -lGLEW -lGL -msse2 -lm -ldl -Wall -pedantic-errors
./bin/debug/hazelholm_c