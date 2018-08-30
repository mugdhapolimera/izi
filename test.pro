fname = 'sortz.txt'
read_data, sortz, fname, COLUMNS = 1, ROWS = 2500

fname = 'sortq.txt'
read_data, sortq, fname, COLUMNS = 1, ROWS = 2500

fname = 'zarr.txt'
read_data, zarr, fname, COLUMNS = 1,ROWS =  50

fname = 'sortpost.txt'
read_data, sortpost, fname, COLUMNS = 1, ROWS = 2500
nz = 50
postz=dblarr(nz)
  for j=0, nz-1 do begin
    postz[j]=int_tabulated(sortq[where(sortz eq zarr[j])], sortpost[where(sortz eq zarr[j])], /SORT)
  endfor

STOP
end