# Procedure for Extracting and Flatfielding IS spectra
#   Use only apnormalized flat
#
# copyright : A.Tajitsu (2013/7/9)
# !!!  It's important to use apall with                      !!!
# !!!          "llimit=(pix) ulimit=(pix+1) ylebel=INDEF"    !!!
# !!!    to extract 1-pixel along a reference aperture.      !!!
#
procedure hdsis_ecf(inimg,outimg)
 string inimg   {prompt= "input image"}
 string outimg  {prompt= "output image\n"}

 bool plot=yes   {prompt= "Plot image and extract manually"}
 int st_x=-12  {prompt ="Start pixel to extract (for plot=no)"}
 int ed_x=12  {prompt ="End pixel to extract (for plot=no)\n"}

 string flatimg {prompt= "ApNormalized flat"}
 string ref_ap  {prompt= "Aperture reference image\n"}

 string badfix="none"  {prompt = 'Fixing method for Bad Pix [none|zero|fixpix]', enum="none|zero|fixpix"}
 real fix_up=0.001 {prompt = 'Upper Limit for Bad Pix in ApNormalized Flat\n'}

struct *list_ecf
begin
#
# variables
#
string inimage, outimage, flat, ref
int low, upp
int i, ysize, i_ord, ans_int
string file_ext, tempspec1, templist
real fmean[100]
string tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
string fix_flt,fix_img,fix_msk,fix_msk1,tmp_flat
int tmp_low, tmp_upp, exptime
real threshold
string ex_flt, ex_img

inimage  = inimg
outimage = outimg

flat = flatimg
ref = ref_ap

ans_int=0

threshold=100000.

#
# start
#

imgets(flat,'i_naxis2')
ysize=int(imgets.value)

imgets(inimage,'EXPTIME')
exptime=int(imgets.value)

if(badfix=="fixpix"){
  fix_flt = mktemp('tmp.hdsis_ecf.fix_flt.')
  fix_msk = mktemp('tmp.hdsis_ecf.fix_msk.')
  fix_msk1 = mktemp('tmp.hdsis_ecf.fix_msk.')
  fix_img = mktemp('tmp.hdsis_ecf.fix_img.')

  imcopy(flat,fix_flt,ver-)
  imcopy(flat,fix_msk1,ver-)
  imcopy(inimage,fix_img,ver-)

  printf("*** Creating Bad Pixel Mask")

  imreplace(fix_msk1, -1, imagina=0.,upper=fix_up,
	    lower=INDEF, radius=0.)
  printf(".")
  imreplace(fix_msk1, 0, imagina=0.,upper=INDEF,
	    lower=fix_up, radius=0.)
  printf(".")
  imarith(fix_msk1,'*',-1,fix_msk)
  printf(".")
  imdelete(fix_msk1)
  printf("done!\n")

  printf("*** Applying Bad Pixel Mask")
  fixpix(fix_flt, fix_msk,linterp=INDEF, cinterp=2, ver-, pixels-)
  printf(".")
  fixpix(fix_img, fix_msk,linterp=INDEF, cinterp=2, ver-, pixels-)
  printf("done!\n")

  ex_img=fix_img
  ex_flt=fix_flt
}
else{
  ex_img=inimage
  ex_flt=flat
}

if(plot){

printf("########## ApEdit for Object frame : %s ##########\n",ex_img)
printf("### Please measure lower & upper for order extraction\n")

apedit(ex_img, referen=ref, aperture="", \
       interac+, find-, recenter-, resize-, edit+, line=INDEF, nsum=10);

list_ecf='database/ap'//ex_img
while(fscan(list_ecf,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)>0)
{
  if(tmp1=='low')
    {
      print(tmp2) | scan(tmp_low);
    }
}
list_ecf='database/ap'//ex_img
while(fscan(list_ecf,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)>0)
{
  if(tmp1=='high')
    {
      print(tmp2) | scan(tmp_upp);
    }
}


printf(">>> Lower Pixel for Extraction (%d) :  ", tmp_low)
if( scan(ans_int) == 0 ){
  print(tmp_low)
  low=tmp_low
}
else{
  print(ans_int)
  low=int(ans_int)
}

printf(">>> Upper Pixel for Extraction (%d) :  ", tmp_upp)
if( scan(ans_int) == 0 ){
  print(tmp_upp)
  upp=tmp_upp
}
else{
  print(ans_int)
  upp=int(ans_int)
}

}
else{
  low=st_x
  upp=ed_x
}

printf("########## Slicing Flat : %s %d -> %d ##########\n",flat, low, upp)

for(i=low;i<upp;i=i+1)
{
     file_ext="_ec"

     if(access(ex_flt//file_ext//i//".fits"))
     {
          imdelete(ex_flt//file_ext//i)
     }

     apall(input=ex_flt,output=ex_flt//file_ext//i,\
      apertur="",\
      format='echelle',reference=ref,profile="",\
      interac-,find-,recente-,resize+,\
      edit-,trace-,fittrac-,extract+,extras-,review-,\
      llimit=i,ulimit=i+1,\
      nsubaps=1, pfit='fit1d', clean-, weights='none',ylevel=INDEF)

     imgets(ex_flt//file_ext//i,'i_naxis2')
     i_ord=int(imgets.value)

     if(badfix=="fixpix"){
       tmp_flat = mktemp('tmp.hdsis_ecf.flat.')

       apall(input=flat,output=tmp_flat,\
        apertur="",\
        format='echelle',reference=ref,profile="",\
        interac-,find-,recente-,resize+,\
        edit-,trace-,fittrac-,extract+,extras-,review-,\
        llimit=i,ulimit=i+1,\
        nsubaps=1, pfit='fit1d', clean-, weights='none',ylevel=INDEF)

       imstat(image=tmp_flat//"[*,2:"//i_ord-1//"]", \
            field='mean', format-,\
            lsigma=3.,usigma=3.,nclip=5,binwidt=0.1) | scan(fmean[i-low+1])

       imdelete(tmp_flat)
     }
     else{
       imstat(image=ex_flt//file_ext//i//"[*,2:"//i_ord-1//"]", \
              field='mean', format-,\
              lsigma=3.,usigma=3.,nclip=5,binwidt=0.1) | scan(fmean[i-low+1])

       if(badfix=="zero"){
           imreplace(ex_flt//file_ext//i, threshold, imagina=0.,upper=fix_up,
  	    lower=INDEF, radius=0.)
       }
     }

     printf(" Pixel = %2d, Mean = %f\n", i, fmean[i-low+1])

}

templist = mktemp('tmp.hdsis_ecf.list.')
printf("########## Slicing Object : %s %d -> %d ##########\n",ex_img, low, upp)
for(i=low;i<upp;i=i+1)
{
     file_ext="_ec"

     printf(" Pixel = %2d\n", i)
     apall(input=ex_img,output=ex_img//file_ext//i,\
      apertur="",\
      format='echelle',reference=ref,profile="",\
      interac-,find-,recente-,resize+,\
      edit-,trace-,fittrac-,extract+,extras-,review-,\
      llimit=i,ulimit=i+1,\
      nsubaps=1, pfit='fit1d', clean-, weights='none',ylevel=INDEF)

     tempspec1 = mktemp('tmp.hdsis_ecf.')
     sarith(ex_img//file_ext//i,"/",ex_flt//file_ext//i,tempspec1)

     file_ext="_ecf"
     sarith(tempspec1,"*",fmean[i-low+1],ex_img//file_ext//i)
     printf("%s\n",ex_img//file_ext//i,>>templist)
     imdelete(tempspec1)
}

printf("########## Combining SubApertures ##########\n")
scombine("@"//templist, outimage, combine="sum", reject="none",\
 group="apertures")
delete(templist)

hedit(outimage,"H_IS_LOW",low,del-,add+,ver-,show-,update+)
hedit(outimage,"H_IS_UPP",upp,del-,add+,ver-,show-,update+)
hedit(outimage,"EXPTIME",exptime,del-,add-,ver-,show-,update+)

if(badfix=="fixpix"){
  templist = mktemp('tmp.hdsis_ecf.list.')
  printf("########## Createing Resultant Mask : Mask_hdsis_%s.fits ###########\n",outimage)
  if(access("Mask_hdsis_"//outimage//".fits")){
     imdelete("Mask_hdsis_"//outimage)
  }
  if(access("Mask_hdsis_"//outimage)){
     delete("Mask_hdsis_"//outimage)
  }
  apall(input=fix_msk,output="Mask_hdsis_"//outimage,\
      apertur="",\
      format='echelle',reference=ref,profile="",\
      interac-,find-,recente-,resize+,\
      edit-,trace-,fittrac-,extract+,extras-,review-,\
      llimit=low,ulimit=upp,\
      nsubaps=1, pfit='fit1d', clean-, weights='none',ylevel=INDEF)
}

for(i=low;i<upp;i=i+1)
{
     file_ext="_ec"
     imdelete(ex_img//file_ext//i)
     imdelete(ex_flt//file_ext//i)
     file_ext="_ecf"
     imdelete(ex_img//file_ext//i)
}

if(badfix=="fixpix"){
  imdelete(fix_flt)
  imdelete(fix_img)
  imdelete(fix_msk)
}

bye
end
