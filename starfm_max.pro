PRO STARFM_max
   
;-----------------------------变量设置----------------------------------------------------------
   
window_size = 25   ;窗口大小，搜索窗口的设置
class_num = 3       ;假定分类数，寻找相同像元时设置
A= 4                ;距离因子,绝对距离转为相对距离设置

;----------------------------------------------------------------------------------------------------------------------   

;----------------------------打开输入数据-------------------------------------------------------------------------------

;open the corse image of the base pair
filename1=dialog_pickfile(title='打开基准modis')
;open the fine image of the base pair
filename2=dialog_pickfile(title='打开基准landsat')
;open the corse image of the pre pair
filename3=dialog_pickfile(title='打开预测modis')
;open the fine image of the pre pair
filename4=dialog_pickfile(title='打开预测landsat')

;-----------------------------读取数据----------------------------------------------------------------
;

;----------------------------打开基准的modis---------------------------------------------------------
envi_open_file,filename1,r_fid=cb_fid
  if cb_fid eq -1 then begin
    envi_batch_exit
    return
  endif
envi_file_query,cb_fid,dims=cb_dims,nl=cb_lines,ns=cb_samples,nb=cb_bands,$
  bnames=cb_bnames,data_type=cb_dt,wl=cb_wl
cb_img=make_array(cb_samples,cb_lines,cb_bands,type=cb_dt)
for nb=0,cb_bands-1 do begin
  cb_img[*,*,nb]=envi_get_data(fid=cb_fid,dims=cb_dims,pos=nb)
endfor  
envi_file_mng,id=cb_fid,/remove

;------------------------打开基准的landsat-----------------------------

envi_open_file,filename2,r_fid=fb_fid
  if fb_fid eq -1 then begin
    envi_batch_exit
    return
  endif
envi_file_query,fb_fid,dims=fb_dims,nl=fb_lines,ns=fb_samples,nb=fb_bands,$
  bnames=fb_bnames,data_type=fb_dt,wl=fb_wl  
fb_img=make_array(fb_samples,fb_lines,fb_bands,type=fb_dt)
for nb=0,fb_bands-1 do begin
  fb_img[*,*,nb]=envi_get_data(fid=fb_fid,dims=fb_dims,pos=nb)
endfor
envi_file_mng,id=fb_fid,/remove

;--------------------------打开预测的modis---------------------------

envi_open_file,filename3,r_fid=cp_fid
if cp_fid eq -1 then begin
  envi_batch_exit
  return
endif
envi_file_query,cp_fid,dims=cp_dims,nl=cp_lines,ns=cp_samples,nb=cp_bands,$
  bnames=cp_bnames,data_type=cp_dt,wl=cp_wl
cp_img=make_array(cp_samples,cp_lines,cp_bands,type=cp_dt)
for nb=0,cp_bands-1 do begin
  cp_img[*,*,nb]=envi_get_data(fid=cp_fid,dims=cp_dims,pos=nb)
endfor
envi_file_mng,id=cp_fid,/remove

;--------------------------------打开预测的landsat----------------------------

envi_open_file,filename4,r_fid=fp_fid
if fp_fid eq -1 then begin
  envi_batch_exit
  return
endif
envi_file_query,fp_fid,dims=fp_dims,nl=fp_lines,ns=fp_samples,nb=fp_bands,$
  bnames=fp_bnames,data_type=fp_dt,wl=fp_wl
fp_img=make_array(fp_samples,fp_lines,fp_bands,type=fp_dt)
for nb=0,fb_bands-1 do begin
  fp_img[*,*,nb]=envi_get_data(fid=fp_fid,dims=fp_dims,pos=nb)
endfor
envi_file_mng,id=fp_fid,/remove

;---------------------------------------------------------------------
;获取投影信息

map_info=envi_get_map_info(fid=fb_fid)
map_info_out=map_info


predicted_data  = fltarr(fp_samples,fp_lines,fp_bands)
;----------------------------------主程序-------------------------------------------------------

t0=systime(1)                  ;初始时间

;------------------------------运行期间参数--------------------------------------
r = (window_size-1)/2 ;移动窗口半径
bianyuan=cp_img+fb_img+cb_img
for  nb=0,fp_bands-1 do begin
  for  sample = r , fp_samples- r - 1 do begin
    for  line = r , fp_lines- r - 1 do begin
   
     
      
      ;获取搜索窗口内的像元
      ;--------------------------------------------------------------------------------------------------------
      data1 = cb_img[(sample-r):(sample+r),(line-r):(line+r),nb]
      data2 = fb_img[(sample-r):(sample+r),(line-r):(line+r),nb]
      data3 = cp_img[(sample-r):(sample+r),(line-r):(line+r),nb]
    
      ;标准差计算，获取窗口内相邻像素
      ;--------------------------------------------------------------------------------------------------------
      similar_pixel_measure = stddev(data2)/class_num
      similar_pixel = intarr(window_size,window_size)
      if similar_pixel_measure ne 0 then begin
        similar_pixel[where(abs(data2-data2[r,r]) lt similar_pixel_measure)] = 1
      endif else begin
        similar_pixel = 1
      endelse
      
      similar_pixel_index = where(similar_pixel eq 1)
      if similar_pixel_index[0] ne -1 then begin
        ;   图像之差
        s_lm = abs(data1[similar_pixel_index] - data2[similar_pixel_index]) > 0.00001
        t_mm = abs(data1[similar_pixel_index] - data3[similar_pixel_index]) > 0.00001
        ;   像素到中心像元的距离 ，并转成相对距离
        sub_index = array_indices([window_size,window_size],similar_pixel_index,/dimensions)
        distance = transpose(1+((sub_index[0,*] - (window_size-1)/2)^2 + (sub_index[1,*] - (window_size-1)/2)^2)^0.5/A)
        ;distance = transpose(1+((sub_index[0,*] - r)^2 + (sub_index[1,*] - r)^2)^0.5*(window_size-1)/2)
        weight = (1/s_lm*t_mm*distance)/total(1/s_lm*t_mm*distance)
        result = total(weight*(data3[similar_pixel_index]+data2[similar_pixel_index]-data1[similar_pixel_index]))
        
      endif else begin
        result = data3[r,r]+data2[r,r]-data1[r,r]
       
      endelse    
        predicted_data[sample,line,nb] = result
      endfor
  endfor
 ; predicted_data[0:r-1,*,nb]=bianyuan[0:r-1,*,nb]
;  predicted_data[fp_samples:r,*,nb]=fb_img[0:r,*,nb]
;  
endfor

print,'ok'



output_filename = 'd:\onedrive\lunwen\test\starfm\2014299.dat'
openw,lun,output_filename,/get_lun
writeu,lun,predicted_data
free_lun,lun
envi_setup_head,fname = output_filename , ns = fp_samples , nl = fp_lines , nb = fp_bands, data_type = fp_dt, offset = 0, interleave = 0, map_info=map_info_out,/write,/open
;plot, fb_img , predicted_data , title='starfm vs. truth' , xtitle='real landsat reflectance' , ytitle='starfm-fused reflectance' , psym = 3 , xrange = [0,1] , yrange = [0,1]
;plot,[0,1],[0,1],xrange = [0,1] , yrange = [0,1],/noerase


print, 'time used:', floor((systime(1)-t0)/3600), 'h',floor(((systime(1)-t0) mod 3600)/60),'m',(systime(1)-t0) mod 60,'s'
END