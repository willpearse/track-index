task :default => [:before, :setup, :gbif, :cru, :calc_index, :figures, :ms, :after]

task :before do
  puts "Running entire analysis from Pearse & Davies (2020)"
  puts "... run `rake --tasks` to see how to run sub-components"
  puts "... open `src/headers.R` to set # of cores for analysis (default 24!)"
end


desc "Install missing R packages"
task :setup do
  `Rscript src/headers.R`
end

desc "Download and process GBIF data"
task :gbif do
  `ruby src/download_gbif.rb`
end

desc "Download and process CRU data"
task :cru do
  `ruby src/download_cru.rb`
  `Rscript src/process-cru.R`
end

desc "Calculating index"
task :calc_index do
  `Rscript src/calculate-index.R`
end

desc "Generating figures"
task :figures do
  `Rscript src/analysis-violin.R`
  `Rscript src/analysis-pca.R`
  `Rscript src/analysis-phylo.R`
  `Rscript src/analysis-traits.R`
end

desc "Building manuscript"
task :ms do 
  `pdflatex ms/ms.tex`
  `pdflatex ms/ms.tex`
  `biber ms/ms`
  `pdflatex ms/ms.tex`
end

task :after do
  puts "... finished!"
end
