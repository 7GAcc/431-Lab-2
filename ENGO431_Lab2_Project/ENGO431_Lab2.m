vHat_27=load("out_vHat_img27.txt");
vHat_28=load("out_vHat_img28.txt");
cert_fid=load("in_calibration_certificate_values2.txt");


xf=cert_fid(2:length(cert_fid),1);
yf=cert_fid(2:length(cert_fid),2);
for i=1:length(vHat_27)/2
   v_xf_27(i)=vHat_27(2*i-1);
   v_yf_27(i)=vHat_27(2*i);
   v_xf_28(i)=vHat_28(2*i-1);
   v_yf_28(i)=vHat_28(2*i);
end


figure(1)
quiver(transpose(xf),transpose(yf),v_xf_27*1000,v_yf_27*1000, 'AutoScale', 'off')
title("Image 27 quiver plot")

figure(2)
quiver(transpose(xf),transpose(yf),v_xf_28*1000,v_yf_28*1000, 'AutoScale', 'off')
title("Image 28 quiver plot")


% ---------Part c---------
test_vHat=load("out_test_vHat.txt");
test_f_coords=load("in_test_reseau.txt");
%test_xvHat=zeros((length(test_vHat/2)));
% test_yvHat=zeros((length(test_vHat/2)));
% test_xf=zeros((length(test_xvHat)));
% test_yf=zeros((length(test_xvHat)));
for i=1:length(test_vHat)/2
    test_xvHat(i)=test_vHat(2*i-1,1);
    test_yvHat(i)=test_vHat(2*i,1);
    
    test_xf(i)=test_f_coords(2*i-1,1);
    test_yf(i)=test_f_coords(2*i,1);
end

figure(3)
quiver(test_xf,test_yf,test_xvHat*1000,test_yvHat*1000, 'AutoScale', 'off')
title("Test quiver plot")

% test_comparator=load("
