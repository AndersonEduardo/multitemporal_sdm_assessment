for /r %i in (*.gz) do gunzip %i

for /r %i in ("landmask") do del %i.asc