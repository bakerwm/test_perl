
## update all packages
for ( lib in .libPaths() ) {
  ## update the packages from all the library paths
  update.packages(lib.loc = lib, ask = FALSE, dependencies = c('Suggests'))
}
