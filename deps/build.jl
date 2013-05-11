using BinDeps

# version of cubature package to use
cubvers="1.0"

tagfile = "installed_vers"
if !isfile(tagfile) || readchomp(tagfile) != cubvers
    if OS_NAME == :Windows
        run(download_cmd("http://ab-initio.mit.edu/cubature/libcubature$WORD_SIZE-$cubvers.dll", "libcubature.dll"))
    elseif OS_NAME == :Darwin
        run(download_cmd("http://ab-initio.mit.edu/cubature/libcubature$WORD_SIZE-$cubvers.dylib", "libcubature.dylib"))
    else
        if !isfile("cubature-$cubvers.tar.gz")
            run(download_cmd("http://ab-initio.mit.edu/cubature/cubature-$cubvers.tgz", "cubature-$cubvers.tar.gz"))
        end
        run(unpack_cmd("cubature-$cubvers.tar.gz", "."))
        cd("cubature-$cubvers") do
            println("Compiling hcubature.c...")
            run(`gcc -fPIC -O3 -c hcubature.c`)
            println("Compiling pcubature.c...")
            run(`gcc -fPIC -O3 -c pcubature.c`)
            println("Linking libcubature...")
            run(`gcc -shared -o ../libcubature.so hcubature.o pcubature.o`)
        end
    end
    run(`echo $cubvers` > tagfile)
end
