#!/usr/bin/ruby

# Change the location (and its match in `headers.R`) to download CRU
# data elsewhere
Dir.chdir "~/bear-share/cru/"

def download_wrapper(cru_id, var, stem=clean=true, stem="cruts.1905011326.v4.03")
  `wget https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.03/#{stem}/#{var}/#{cru_id}`
  `gunzip #{cru_id}`
  if clean then `rm #{cru-id}` end
  Dir.chdir '..'
end

download_wrapper("cru_ts4.03.1901.2018.cld.dat.nc.gz", "cld")
download_wrapper("cru_ts4.03.1901.2018.dtr.dat.nc.gz", "dtr")
download_wrapper("cru_ts4.03.1901.2018.frs.dat.nc.gz", "frs")
download_wrapper("cru_ts4.03.1901.2018.pet.dat.nc.gz", "pet")
download_wrapper("cru_ts4.03.1901.2018.pre.dat.nc.gz", "pre")
download_wrapper("cru_ts4.03.1901.2018.tmn.dat.nc.gz", "tmn")
download_wrapper("cru_ts4.03.1901.2018.tmp.dat.nc.gz", "tmp")
download_wrapper("cru_ts4.03.1901.2018.tmx.dat.nc.gz", "tmx")
download_wrapper("cru_ts4.03.1901.2018.vap.dat.nc.gz", "vap")
download_wrapper("cru_ts4.03.1901.2018.wet.dat.nc.gz", "wet")
