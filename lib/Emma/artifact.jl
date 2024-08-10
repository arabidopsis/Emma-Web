using Tar, Inflate, SHA
name = "Emma2-model"
url = "https://github.com/ian-small/emma-models/archive/refs/tags/v1.0.0.tar.gz"
filename = basename(url)
rm(filename, force=true)

run(pipeline(`wget $(url)`, stderr=devnull))
println("[$(name)]")
println("git-tree-sha1 = \"$(Tar.tree_hash(IOBuffer(inflate_gzip(filename))))\"")
println("")
println("    [[$(name).download]]")
println("    url = \"$(url)\"")
println("    sha256 = \"$(bytes2hex(open(sha256, filename)))\"")
rm(filename, force=true)
