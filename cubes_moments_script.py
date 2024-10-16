import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import SigmaClip
from photutils import MADStdBackgroundRMS
import numpy as np
#example; from now on I am going to consider i=0
gal=['ESO320-G030', 'IC4518E']



#To identify the CO(2–1) emission in each channel of the data cube, we selected pixels with fluxes > 5sigma.
# So we clip every channel to create a clipped data cube
for i in range(len(gal)):
    i=0
    new_cube=[]
    sig_ch=[]
    hdu=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'/CO21_image_wtPB.fits')[0]   # datacube without primary beam correction
    x=np.shape(hdu) #  print(x)  (1, 325, 450, 450) 325 channels
    header=hdu.header
    xx=x[1] # number of channels
    for j in range(xx):  # for each channel I want to estimate the sigma to make a clipping
        channel_data = hdu.data[0,j,:,:]
        oo=channel_data[np.isnan(channel_data)==False] # remove nan values
        sigma_clip = SigmaClip(sigma=3.0) # 3 is the default value for sigma clipping.
        bkgrms = MADStdBackgroundRMS(sigma_clip) # to calculate the background RMS in an array as using the median absolute deviation (MAD).
        bkgrms_value_channel = bkgrms.calc_background_rms(oo)
        sigma=bkgrms_value_channel
        print(sigma, j)
        sig_ch.append(sigma) # to create an array with the values of sigma for each channel of the cube
    
        channel_data[channel_data<5*sigma]=0 # I set to 0 all values minor than the value of 5*sigma, I choose 5*sigma because a lower value of sigma introduces some noise pixels
        new_cube.append(channel_data)
 
    
    data=new_cube
    hdu=fits.PrimaryHDU(data=new_cube, header=header) # new data cube
    HST_flux='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_cube-wtPB-mod_5sig.fits' # name of the new datacube 'mod=modified' and '5sig=5 sigma' 

    fits.writeto(HST_flux,data,header, overwrite=True)


#### To obtain the clipped cube with the Primary beam correction
for i in range(len(gal)):

    hdu1=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_cube-wtPB-mod_5sig.fits')[0]
    x1=np.shape(hdu1.data)
    hdu1.data[hdu1.data==0]=np.nan

    hdu2=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/CO21.flux.fits')[0] # Primary Beam image

    x2=np.shape(hdu2.data)

    hf=np.divide(hdu1.data,hdu2.data) 
    hdu=fits.PrimaryHDU(data=hf, header=hdu1.header)

    HST_flux='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_cube-PBcor-mod_5sig.fits'
 
    fits.writeto(HST_flux,hf,header=hdu1.header, overwrite=True)



## To avoid including noise spikes, we remove the isolated channels with emission 
for i in range(len(gal)):
    hdul=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_cube-PBcor-mod_5sig.fits')[0]


    dat=np.squeeze(hdul.data)
    dat[np.isnan(dat)]=0
    x1=np.shape(dat)
    x10=x1[0]
    x11=x1[1]
    x12=x1[2]

    b=[]
    for j in range(x11):   
        for k in range(x12):
            for m in range(x10):
                aa=dat[m,j,k] 
                b.append(aa)
            for l in range(1, x10-1, 1):
                if b[l-1]==0 and b[l+1]==0:
                    b[l]=0.0
                else:
                    b[l]=b[l]
            dat[:,j,k]=b
            print(j, k)
            b=[] 
    data=dat
    data[data==0]=np.nan
    hdu=fits.PrimaryHDU(data=data, header=hdul.header)
    HST_flux='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_cube-PB-mod_5sig_no_noise.fits'

    fits.writeto(HST_flux,data=data,header=hdul.header, overwrite=True)



# To create moment 0 from the definition. This moment 0 comes from the previous steps (a clipped data cube + no isolated channels with emission)
#dv is the spectral resolution
dv=np.array([10.159116, 10.159116]) # example # dv is the channel width
for i in range(len(gal)):    

    hdul=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_cube-PB-mod_5sig_no_noise.fits')[0]
    hdul.data[np.isnan(hdul.data)] = 0
    y=np.squeeze(hdul.data) # to  check what is the meaning of this
    x=np.shape(y)
    x0=x[0]
    x1=x[1]
    x2=x[2]
    v=np.zeros((x1,x2))
    for m in range(x0):
        v=np.add(v,y[m,:,:])
    v=v*dv[i]
    v[v==0]=np.nan

    fig, ax = plt.subplots()
    plt.imshow(v,origin='lower')
    plt.savefig('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_mom0-PBcor-mod_5sig_no_noise.pdf')
 
    hdu=fits.PrimaryHDU(data=v, header=hdul.header)
    HST_flux='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_mom0-PBcor-mod_5sig_no_noise.fits'

    fits.writeto(HST_flux,data=v,header=hdul.header, overwrite=True)
    


# now I want to create the map with the number of channels with emission that each pixel has
for i in range(len(gal)):
    position_x=[]
    position_y=[]
    value_channel=[]
   
    hdul=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_cube-PB-mod_5sig_no_noise.fits')[0]
    hdul.data[np.isnan(hdul.data)] = 0
    data=np.squeeze(hdul.data)
    x=np.shape(data)
    x0=x[0]
    x1=x[1]
    x2=x[2]
    w=np.zeros((x1,x2))  
    for j in range(x1):
        for k in range(x2):
            for l in range(x0):
                if data[l,j,k] > 0:
                    w[j,k]=w[j,k]+1
                    # print(j,k)
    w[w==0]=np.nan
    maxi=np.nanmax(w)
    mini=np.nanmin(w)
    print(gal,maxi,mini)
    fig, ax = plt.subplots()
    plt.imshow(w,origin='lower',cmap=plt.cm.get_cmap('jet', maxi-1)) 
    plt.colorbar()
    plt.clim(mini, maxi);
    plt.savefig('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'-channels_pbcor_created_5sig.pdf') # this is the figure of the map with the number of channels per pixel 
 
    w[w<=2]=np.nan
    fits.PrimaryHDU(data=w, header=hdul.header)
    HST_flux='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_number_of_channels_in_moment0_5sig-2pix.fits'# name of the fits
 
    fits.writeto(HST_flux,w,hdul.header, overwrite=True)

    mx=np.nanmax(w)
    fig, ax = plt.subplots()
    plt.imshow(w,origin='lower',cmap=plt.cm.get_cmap('jet', mx-2)) # I have to check this
    plt.colorbar()
    plt.clim(2, mx);
   
    plt.savefig('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'-number_of_channels_pbcor_created_5sig_pix2.pdf') # no se que hago aqui

    w[np.isnan(w)]=0
    w[w>2]=1
    fig, ax = plt.subplots()
    plt.imshow(w,origin='lower')
    plt.colorbar()
    plt.savefig('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'-imagen_binary-5sig.pdf') # the binary image

    fits.PrimaryHDU(data=w, header=hdul.header)
    HST_flux='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_binary-5sig-2pix.fits' # binary fits

    fits.writeto(HST_flux,w,hdul.header, overwrite=True)


    hdu2=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_mom0-PBcor-mod_5sig_no_noise.fits')[0] # import the moment 0
                                                                        
    data2=np.squeeze(hdu2.data)
    fig, ax = plt.subplots()
    plt.imshow(data2,origin='lower')
    plt.colorbar()


    dat=data2/w # I devide the moment 0 by the binary to obtain a moment 0 where the pixels are obtained from more than 2 channels of emission
    dat[dat==np.inf]=np.nan
    fig, ax = plt.subplots()
    plt.imshow(dat,origin='lower')
    plt.colorbar()
    plt.savefig('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'-5sigma_2pix.pdf') # This is the plot of the final moment 0

    fits.PrimaryHDU(data=dat, header=hdul.header)
    HST_flux='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_mom0-PBcor-mod_5sig-2pix_no_noise.fits' # This is the final moment 0

    fits.writeto(HST_flux,dat,hdul.header, overwrite=True)


# I add a channel in the beginning and in the end of each pixel with emission
for i in range(len(gal)):

    hdul5=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_cube-PB-mod_5sig_no_noise.fits')[0]

    dat=np.squeeze(hdul5.data)
    data_new=dat
    dat[np.isnan(dat)]=0
    dat[dat==np.inf]=0
    hdul=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/CO21.image.pbcor.fits')[0]

    dat1=np.squeeze(hdul.data)
    dat1[np.isnan(dat1)]=0
    dat1[dat1==np.inf]=0

    x1=np.shape(dat)
    x10=x1[0]
    x11=x1[1]
    x12=x1[2]
    b=[]
    b_new=np.zeros(x10)
    b1=[]
    for j in range(x11):   
        for k in range(x12):
            for m in range(x10):
                aa=dat[m,j,k] 
                b.append(aa)
                aa1=dat1[m,j,k] 
                b1.append(aa1)
            for l in range(1, x10-1, 1):
                print(l)
                if b[l]==0 and (b[l-1] != 0 or b[l+1] != 0):
                    b_new[l]=b1[l]
                    print(b_new[l])
                else:
                    b_new[l]=b[l]
            data_new[:,j,k]=b_new
            print(j, k)  
            b=[]
            b1=[]
            b_new=np.zeros(x10)
    data_new[data_new==0]=np.nan
    fits.PrimaryHDU(data=data_new, header=hdul5.header)
    HST_flux='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_cube-PB-mod_5sig_no_noise_plus-channels.fits'
    fits.writeto(HST_flux,data=data_new,header=hdul5.header, overwrite=True)


# Moment 0 with the previous conditions: 
for i in range(len(gal)):  
    i=0
    hdul=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_cube-PB-mod_5sig_no_noise_plus-channels.fits')[0]

    hdul.data[np.isnan(hdul.data)] = 0
    x=np.shape(hdul.data)
    x0=x[0]
    x1=x[1]
    x2=x[2]
    v=np.zeros((x1,x2))

    for m in range(x0):
        v=np.add(v,hdul.data[m,:,:])
    v=v*dv[i]
    v[v==0]=np.nan

    fig, ax = plt.subplots()
    plt.imshow(v,origin='lower')
    plt.savefig('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_mom0-PBcor-mod_5sig-no_noise_plus-channels.pdf') 

    fits.PrimaryHDU(data=v, header=hdul.header)
    HST_flux='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_mom0-PBcor-mod_5sig-no_noise_plus-channels.fits'

    fits.writeto(HST_flux,data=v,header=hdul.header, overwrite=True)




"""
moment -1
"""
from astropy import wcs
from astropy import units as u
for i in range(len(gal)):
    hdul=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/CO21.image.pbcor.fits')[0]      

    headeri=hdul.header
    bmaj=headeri['BMAJ']
    bmin=headeri['BMIN']
    pa=headeri['BPA']
    hdu=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/CO21_image_wtPB.fits')[0] 

    headerr=hdu.header
    data=np.squeeze(hdul.data) 
    data[np.isnan(data)] = 0
    data_final=np.amax(data,axis=0)   
    data_final[data_final==0]=np.nan
    headerr['BUNIT']= "Jy/beam"
    headerr['BMAJ']=bmaj    
    headerr['BMIN']=bmin
    headerr['BPA']=pa 
    hdu=fits.PrimaryHDU(data=data_final, header=headerr)
    HST_flux='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_mom1n-PBcor.fits'

    fits.writeto(HST_flux,data_final,headeri, output_verify='ignore', overwrite=True)
    hdu2=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/CO21.flux.fits')[0]

    dataa=np.squeeze(hdu2.data)

    hf=np.divide(data_final,dataa[1,:,:])
    hdu=fits.PrimaryHDU(data=hf, header=headerr)
    HST_flux='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_mom1n-PBcor.fits'

    fits.writeto(HST_flux,data_final,headeri, output_verify='ignore', overwrite=True)




"""""

Rest of the moments: moment 1  

"""""


for i in range(len(gal)):
    alma='/Volumes/Maria/cubes/0.data/'+gal[i]+'/CO21.image.pbcor.fits'      

    header = fits.getheader(alma)  

    hdulist=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_cube-PB-mod_5sig_no_noise_plus-channels.fits')[0]
  
    head=hdulist.header
    data=hdulist.data
    data=np.squeeze(data)
    data[np.isnan(data)] = 0
    w = wcs.WCS(hdulist.header)
    freq_wcs = w.sub([3]) # this gives information about the third dimension
    freq_axis = freq_wcs.all_pix2world(np.arange(data.shape[0]), 0)[0]*u.Hz
  #  freq_to_vel = u.doppler_relativistic(freq_wcs.wcs.restfrq*u.Hz)
    freq_to_vel = u.doppler_radio(freq_wcs.wcs.restfrq*u.Hz)
    vel_axis = freq_axis.to(u.km / u.s, equivalencies=freq_to_vel)
    # here I obtain the velocity
    x=np.shape(data)
    x0=x[0]
    x1=x[1]
    x2=x[2]     
    num=np.zeros((x1,x2))
    den=np.zeros((x1,x2))
    new_cube=[]
    for n in range(x0): # I create the new cube
        channel_data=data[n,:,:]
        channel_data=channel_data*vel_axis[n] # I put to 0 all values minor than the condition
        channel_data[np.isnan(channel_data)] = 0
        new_cube.append(channel_data)
    cube_new=np.array(new_cube)
    for n in range(x0):
        num=np.add(num,cube_new[n,:,:])
        den=np.add(den,data[n,:,:])
    num[num==0]=np.nan
    den[den==0]=np.nan
    mom1=num/den
    fig, ax = plt.subplots()
    plt.imshow(mom1,origin='lower', cmap=plt.cm.jet, vmin=np.nanmin(mom1), vmax=np.nanmax(mom1))
    plt.colorbar()
    plt.savefig('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'-mom1_vradio_created_5sig-no_noise-plus-channels.pdf')

    head['BUNIT']= "km/s"
    head['BMAJ']=bmaj
    head['BMIN']=bmin
    hdu=fits.PrimaryHDU(data=mom1, header=head)
    HST_fluxx='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'-mom1_vradio-PBcor-5sig-no_noise-plus-channels.fits'

    fits.writeto(HST_fluxx,data=mom1,header=head, overwrite=True)
    # If I don't want to consider the pixels with less than 3 channels with emission
    binario=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_binario-5sig-2pix.fits')[0]
    
    data_binario=binario.data

    mom1_final=mom1/data_binario
    mom1_final[mom1_final==0]=np.nan
    mom1_final[mom1_final==np.inf]=np.nan
    fig, ax = plt.subplots()
    plt.imshow(mom1_final,origin='lower', cmap=plt.cm.jet, vmin=np.nanmin(mom1_final), vmax=np.nanmax(mom1_final))
    plt.colorbar()
    plt.savefig('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'-mom1_vradio_created_5sig-2pix-no_noise-plus-channels.pdf')

    hdu=fits.PrimaryHDU(data=mom1_final, header=head)
    HST_flux='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'-mom1_vradio-PBcor-5sig-2pix-no_noise-plus-channels.fits'

    fits.writeto(HST_flux,data=mom1_final,header=head, overwrite=True)

      

"""""

Rest of the moments: moment 2

"""""


for i in range(len(gal)):
    alma='/Volumes/Maria/cubes/0.data/'+gal[i]+'/CO21.image.pbcor.fits'

    header = fits.getheader(alma)  
    hdulist=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_cube-PB-mod_5sig_no_noise_plus-channels.fits')[0]

    head=hdulist.header
    data=hdulist.data
    data=np.squeeze(data)
    data[np.isnan(data)] = 0
  #  print(np.shape(data))
    w = wcs.WCS(hdulist.header)
    freq_wcs = w.sub([3]) # this gives information about the third dimension
#    lambs = head['CRVAL3'] + head['CDELT3']*(np.linspace(1, head['NAXIS3'], head['NAXIS3']) - head['CRPIX3'])
    # lambs is in Hz
    freq_axis = freq_wcs.all_pix2world(np.arange(data.shape[0]), 0)[0]*u.Hz
#    freq_to_vel = u.doppler_relativistic(freq_wcs.wcs.restfrq*u.Hz)
    freq_to_vel = u.doppler_radio(freq_wcs.wcs.restfrq*u.Hz)
    vel_axis = freq_axis.to(u.km / u.s, equivalencies=freq_to_vel)
    # here I obtain the velocity
    x=np.shape(data)
    x0=x[0]
    x1=x[1]
    x2=x[2]     
    num=np.zeros((x1,x2))
    den=np.zeros((x1,x2))
    new_cube=[]
    moment2 = fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'-mom1_vradio-PBcor-5sig-2pix-no_noise-plus-channels.fits')[0]

    head2=moment2.header
    data2=moment2.data
    vel=np.array(vel_axis)
    for n in range(x0): # I create the new cube
        channel_data=data[n,:,:]
        channel_data=channel_data*(vel[n]-data2)**2 # I put to 0 all values minor than the condition
        channel_data[np.isnan(channel_data)] = 0
        new_cube.append(channel_data)
    cube_new=np.array(new_cube)
    for n in range(x0):
        num=np.add(num,cube_new[n,:,:])
        den=np.add(den,data[n,:,:])
    num[num==0]=np.nan
    den[den==0]=np.nan
    abcd=num/den
    mom2=np.sqrt(num/den)
    fig, ax = plt.subplots()
    plt.imshow(mom2,origin='lower', cmap=plt.cm.jet, vmin=np.nanmin(mom2), vmax=np.nanmax(mom2))
 #   plt.title(my_galaxies[i])
    plt.colorbar()
    plt.savefig('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'-mom2_vradio_created_5sig-no_noise-plus-channels.pdf')

    head['BUNIT']= "km/s"
    head['BMAJ']=bmaj
    head['BMIN']=bmin
    head['BPA']=pa

    hdu=fits.PrimaryHDU(data=mom2, header=head)
    HST_fluxx='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'-mom2_vradio-PBcor-5sig-no_noise-plus-channels.fits'


    fits.writeto(HST_fluxx,data=mom2,header=head, overwrite=True)
    binario=fits.open('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'_binario-5sig-2pix.fits')[0]

    data_binario=binario.data

    mom2_final=mom2/data_binario
    mom2_final[mom2_final==0]=np.nan
    mom2_final[mom2_final==np.inf]=np.nan
    fig, ax = plt.subplots()
    plt.imshow(mom2_final,origin='lower', cmap=plt.cm.jet, vmin=np.nanmin(mom2_final), vmax=np.nanmax(mom2_final))
  #  plt.title(my_galaxies[i])
    plt.colorbar()   
    plt.savefig('/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'-mom2_vradio_created_5sig-2pix-no_noise-plus-channels.pdf')
  
    head['BUNIT']= "km/s"
    head['BMAJ']=bmaj
    head['BMIN']=bmin
    head['BPA']= pa

    hdu=fits.PrimaryHDU(data=mom2_final, header=head)
    HST_flux='/Volumes/Maria/cubes/0.data/'+gal[i]+'/'+gal[i]+'-mom2_vradio-PBcor-5sig-2pix-no_noise-plus-channels.fits'

    fits.writeto(HST_flux,data=mom2_final,header=head, overwrite=True)



