#
# produce input files for splitting
#
#
# run like
#
#
# echo $depths $densities | gawk -v tensor_dir=tensors/ -v density_bot=3.513 \
#                                -v vp_bot=8.809 -v vs_bot=4.734 -v rayp=0.01834 \
#                                -v az=90 -f produce_splitting_layers.awk 
#
# Modified after a script by Thorsten Becker
#
BEGIN{
  # variables 
  if(tensor_dir == "")
    print("error, tensor_dir needs to be defined") > "/dev/stderr";
  if(az == "")
    print("error, az needs to be defined") > "/dev/stderr";
  if(vp_bot == "")  
    print("error, vp_bot needs to be defined") > "/dev/stderr";
  if(vs_bot == "")
    print("error, vs_bot needs to be defined") > "/dev/stderr";
  if(density_bot == "")
    print("error, density_bot needs to be defined") > "/dev/stderr";
  if(rayp == "")
    print("error, rayp needs to be defined") > "/dev/stderr";
  if(idisp == "")
    idisp = 1;
  if(outfile == "")
    outfile = "stdout";
}
{
  if((NF-1)%2 != 0){
    print("error, paste depths and densities",NF%2) > "/dev/stderr";
  }else{
    nla = (NF-1)/2;
for(i=1;i<=nla;i++){
      d[i] = $i;		# depths
      r[i] = $(i+nla);		# density
      #print(d[i],r[i]) > "/dev/stderr";
    }
  shift=$(1+2*nla)
  #print(shift) > "/dev/stderr";
  #d[nla+1]=2*(100-shift)-d[nla];
  d[nla+1]=2*(bottom-shift)-d[nla];
  }
}
END{
  tnla = nla + 2;
  ztotal = 0;
  printf("%i\t\t\t!number of layers (normally 9)\n",tnla);
  printf("0\t\t\t!-------------- air layer to model free surface\n");
  printf("8            		!menu, 8=simple % anisotropy, non-poisson solid\n");
  printf("0.0001 0.0001 0.0000001 !Vp, Vs, % Vp anisotropy\n");
  printf("0            		!4-theta factor\n");
  printf("2            		!symmetry axis, 1=slow, 2=fast\n");
  printf("0.00001        		!density\n");
  printf("1 90         		!rotation axis, angle \n");
  printf("3 0          		!rotation axis, angle\n");
  for(i=1;i<=nla;i++){
    if(i == 1){
      lthick = (d[i+1] + d[i])/2; # layer thickness
    }else if(i >1)#!= nla)
      lthick = (d[i+1] + d[i])/2 - (d[i] + d[i-1])/2;
    #if(i != 1)
      zlevel=d[i];		# bottom level of layer
     # zlevel += lthick;
    ztotal += lthick;
    printf("%g\t\t\t!--- layer thickness - tensor %i at %g km (%g - %g)\n",lthick,i,zlevel,
	   ztotal-lthick,ztotal);
    if(az==0)
	printf("%g\t\t\t!--- layer thickness - tensor %i at %g km (%g - %g)\n",lthick,i,zlevel,
	   ztotal-lthick,ztotal) > "/dev/stderr";
    if(zlevel != d[i])
	print("error with depths",zlevel,d[i]) > "/dev/stderr";
    printf("9		        !menu, 9 = read full 81 component tensor from file\n");
    printf("%s/depth_%g.cijkl\n",tensor_dir,d[i]+shift);
    printf("1.		        !scale factor for elastic coeffs; need GPa/rho[g/cm^3]\n");
    printf("%g\t\t\t!density\n",r[i]);
    printf("1 90        !rotation axis, angle (change anisotropy tilt here)\n");
    printf("3 %g 	        !rotation axis, angle (change incidence azimuth here)\n",az);
  }
  laststep = 20;#lthick/2;
  ztotal += laststep;
  printf("%g         		!---- bottom layer for incident polarization at %g km\n",
	 laststep,ztotal);
  if(az==0)
     printf("%g         		!---- bottom layer for incident polarization at %g km\n",
	 laststep,ztotal) > "/dev/stderr";
  printf("8            		!menu, 8=simple % anisotropy, non-poisson solid\n");
  printf("%g %g 0.000001\t\t	!Vp, Vs, % Vp anisotropy\n",vp_bot,vs_bot);
  printf("0            		!4-theta factor\n");
  printf("2            		!symmetry axis, 1=slow, 2=fast\n");
  printf("%g\t\t\t!density\n",density_bot);
  printf("1 90        		!rotation axis, angle (change anisotropy tilt here)\n");
  printf("3 0         		!rotation axis, angle (change incidence azimuth here, -999 will loop)\n");
  printf("%g\t\t\t! ray parameter\n",rayp);
  printf("0.0 25  .006103515625   !frequency min,max,spacing\n");
  printf("1            		!grt/c matrix, 1=tran_u, 2=ref_d\n");
  printf("%s\n",outfile);
  #if(az==0) printf("%s\n",outfile) > "/dev/stderr";
  printf("1            		!number of output depths\n");
  printf("1            		!layer number for output\n");
  printf("%i         		!1=short display, 2=long display 3: seis 4,5,6,7: P,SKS,SI,SK2, +10, +20 for splitting\n",idisp);
}
