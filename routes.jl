using Genie, Genie.Router, Genie.Renderer.Html, Genie.Requests
using FASTX
using BioSequences
using ZipFile
include("lib/Emma/src/Emma.jl")


route("/") do
  serve_static_file("annotate.html")

end

route("/", method = POST) do

  #remove existing zip files and directories
  filelist = ["public/gffdir","public/gbdir", "public/gbdir","public/fadir","public/gb.zip","public/gff.zip","public/fa.zip", "genomes"]
  for file in filelist
    isfile(file) && rm(file)
    isdir(file) && rm(file, recursive = true)
  end

  #Html for the header
  header = """<head>
  <title>Emma</title>
  <link rel="stylesheet" type="text/css" href="css/genie/style.css">
  </head>
  <body>
  <div class="green-banner">
    <h1 class="centered">Emma</h1>
    <h3 class="centered">Vertebrate mitogenome annotator</h3>
  </div>
  <div class="tab-container">
    <div class="first-tab">Emma</div>
    <img class="tab-icon"src="img/pen.png">
    <a class="tab" href="annotate.html">Annotate</a>
    <img class="tab-icon"src="img/book.png">
    <a class="tab" href="info.html">Info</a>
    <img class="tab-icon"src="img/phone.png">
    <a class="tab" id="contact.html">Contact</a>
  </div>"""
  
  if infilespayload(:yourfile)
    write(filespayload(:yourfile))
    file = filespayload(:yourfile).data
    if haskey(postpayload(), :myCheckbox)
      fastaout = "public/result.fa";
    else
      fastaout = nothing
    end

    #get file from app
    temp_path = tempname();
    open(temp_path, "w") do f
      write(f, filespayload(:yourfile).data)
    end


    println(nthreads())
    try #First, check the input to see if it is a zip
      mkdir("genomes")
      global ARGS = ["--gff", "public/gffdir", "--gb", "public/gbdir"]
      if fastaout != nothing
        push!(ARGS, "--fa", "public/fadir")
      end
      push!(ARGS, "genomes")
      r = ZipFile.Reader(temp_path)
      for f in r.files
        write("genomes/$(f.name)", f)
      end
      close(r)
      mkdir("public/gffdir")
      mkdir("public/gbdir")
      mkdir("public/fadir")
      main()

      gff = ZipFile.Writer("public/gff.zip")
      for f in readdir("public/gffdir")
        w = ZipFile.addfile(gff, f, method=ZipFile.Deflate);
        f = join("public/gffdir/$f")
        write(w, read(f, String))
      end
      close(gff)
      gb = ZipFile.Writer("public/gb.zip")
      for f in readdir("public/gbdir")
        w = ZipFile.addfile(gb, f, method=ZipFile.Deflate);
        f = join("public/gbdir/$f")
        write(w, read(f, String))
      end
      close(gb)
      gff = ZipFile.Writer("public/gff.zip")
      for f in readdir("public/gffdir")
        w = ZipFile.addfile(gff, f, method=ZipFile.Deflate);
        f = join("public/gffdir/$f")
        write(w, read(f, String))
      end
      close(gff)
      gb = ZipFile.Writer("public/fa.zip")
      for f in readdir("public/fadir")
        w = ZipFile.addfile(gb, f, method=ZipFile.Deflate);
        f = join("public/fadir/$f")
        write(w, read(f, String))
      end
      close(gb)

      html(header *
      """<div class="card-container">
        <div class="card">
          <a class="card-header" href=gff.zip download=gff.zip>Download gff</a>
          <a class="card-header" href=gb.zip download=gb.zip>Download gb</a>
          <a class="card-header" href=fa.zip download=fa.zip>Download fasta</a>
        </div>
      </div>""")
    catch ex #If something fails, expect single fasta file
      global ARGS = ["--gff", "public/result.gff", "--gb", "public/result.gb", "--svg", "public/result.svg"]
      if fastaout != nothing
        push!(ARGS, "--fa", "public/result.fa")
      end
      gffdownload = "result.gff"
      gbdownload = "result.gb"
      push!(ARGS, temp_path)
      main()

      html(header *
      """<div class="card-container">
        <div class="card">
          <a class="card-header" href=result.gff download=gff.zip>Download gff</a>
          <a class="card-header" href=result.gb download=gb.zip>Download gb</a>
          <a class="card-header" href=result.fa download=test.fa>Download fasta</a>
        </div>
        </div>
        <div class="card-container">
          <img class="card" src="result.svg">
        </div>
        """)
    end
  else
    println("not a file")
  end
end