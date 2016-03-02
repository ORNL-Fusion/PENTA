pro plot_penta_results

fpath='/home/jjl/PENTA/Trunk/TESTS/mytest/'
fname = fpath+'fluxes_vs_Er'
;openr,lun,fname,/get_lun
;free_lun,lun

data = read_ascii(fname,data_start=2)
data = data.field1

window,/free
c=mycolors()
plot,data[1,*],data[2,*],xrange=[220.,240.],xstyle=1
oplot,data[1,*],data[3,*],color=c.blue



end
