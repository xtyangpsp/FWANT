BEGIN {
    ROOT="/depot/xtyang/data/projects/noahtrupin/LSQR_test/example"
}
{ gsub(/\.\.\//, ROOT "/"); print }
