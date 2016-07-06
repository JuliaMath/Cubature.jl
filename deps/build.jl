using BinDeps
using Compat
import Base.Sys.WORD_SIZE

# version of cubature package to use
cubvers="1.0.2"

tagfile = "installed_vers"
if !isfile(tagfile) || readchomp(tagfile) != "$cubvers $WORD_SIZE"
    info("Installing Cubature $cubvers library...")
    if is_windows()
        run(download_cmd("http://ab-initio.mit.edu/cubature/libcubature$WORD_SIZE-$cubvers.dll", "libcubature.dll"))
    elseif is_apple()
        run(download_cmd("http://ab-initio.mit.edu/cubature/libcubature$WORD_SIZE-$cubvers.dylib", "libcubature.dylib"))
    else
        if !isfile("cubature-$cubvers.tar.gz")
            run(download_cmd("http://ab-initio.mit.edu/cubature/cubature-$cubvers.tgz", "cubature-$cubvers.tar.gz"))
        end
        run(unpack_cmd("cubature-$cubvers.tar.gz", ".", ".gz", ".tar"))
        cd("cubature-$cubvers") do
            println("Compiling hcubature.c...")
            run(`gcc -fPIC -O3 -c hcubature.c`)
            println("Compiling pcubature.c...")
            run(`gcc -fPIC -O3 -c pcubature.c`)
            println("Linking libcubature...")
            run(`gcc -shared -o ../libcubature.so hcubature.o pcubature.o`)
        end
    end
    open(tagfile, "w") do f
        println(f, "$cubvers $WORD_SIZE")
    end
else
    info("Cubature $cubvers is already installed.")
end
