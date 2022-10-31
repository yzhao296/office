%intergral of z at kz = 0--------------------------------------
clear all
clc
close all
spec = load('spec.mat');
spec = cell2mat(struct2cell(spec));

path = 'C:\Users\User\Desktop\yunlei\ecoli_bin2\';
list=dir([path, '*.tif']);
% temp = imread(dc(1).name);
% tempsize = size(temp);
% kai = zeros(:,:,length(dc));
for i = 1 :1: length(list) % loop through all your images to resize

    kai(:,:,i) = imread([path,list(i).name]);
% newimage=imresize(image, 2);

end
% kaiq

% % [Xq,Yq,Zq] = meshgrid(.1:.25:10,-3:.25:3,-3:.25:3);
% kai = interp3(kai,'cubic',10);
[Li,Lj,Lk]=size(kai);
% kai = 4*kai;

refractive_index = 1.3;
theta = 7*pi/18;
phi_scan_number = 36;
sa_index = 4;



%under 10X esitimation, unit m---------------------
kai_size1 = 6e-8;
kai_size2 = 5e-8;
%particle setting---------------------------------
particle_size1 = Li;
particle_size2 = Lj;
particle_size3 = Lk;

%image plane setting-------------------------
image_size = 201;
image_center_index = round(image_size/2);


%scattering potential kai------------------
% kai = zeros(particle_size1,particle_size2,particle_size2);
particle_center_index1 = round(particle_size1/2);
particle_center_index2 = round(particle_size2/2);
particle_center_index3 = round(particle_size3/2);
z_plus = 0;
z_minus = 0;



ffted_kai = fftshift(fftn(ifftshift(kai)));

figure(4)
imagesc(kai(:,:,particle_center_index3));
set(gca,'dataAspectRatio',[1,1,5/6])
set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
set(gca,'visible','off')
figure(5)
isosurface(abs(kai),1e-8);
set(gca,'dataAspectRatio',[1,1,6/5])
set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
set(gca,'visible','off')

figure(6)
imagesc(log(abs(ffted_kai(:,:,51))));
set(gca,'dataAspectRatio',[1,1,1])
figure(7)
isosurface(log(abs(ffted_kai)),1e-4);
set(gca,'dataAspectRatio',[1,1,1])
figure(8)
isosurface(abs(ffted_kai),1e1);
set(gca,'dataAspectRatio',[1,1,1])
temp=fftshift(ifftn(ifftshift(ffted_kai)));
figure(9)
isosurface(abs(temp),1e-4);
set(gca,'dataAspectRatio',[1,1,1])




Aper = zeros(image_size,image_size);

for i = 1:image_size
    for j = 1:image_size
        for k = 1:image_size
             if ((i-image_center_index))^2 + ((j-image_center_index))^2 < (0.6*(k-image_center_index))^2
                 Aper(i,j,k)=1;
             end
        end
    end
end
% figure(10)
% isosurface(Aper,1e-4);
% set(gca,'dataAspectRatio',[1,1,1])


% Aper = ones(image_size,image_size);


%alpha less than 1e12 contributes almost no influence because kai is binary
%0 or 1
alpha_f = 10;
alpha_b = 10;



fI_plus = zeros(image_size,image_size);
fI_minus = zeros(image_size,image_size);
fI_f = zeros(image_size,image_size);
fI_b = zeros(image_size,image_size);
fI_plus_sa = zeros(image_size,image_size);
fI_minus_sa = zeros(image_size,image_size);
fI_f_sa = zeros(image_size,image_size);
fI_b_sa = zeros(image_size,image_size);


for n = 160:50:1760
    tic
    %parameters---------------------------------
%     lamda = 1e-9*spec(752,1);
    lamda = 1e-9*spec(n,1);
    k_max = 8*pi*refractive_index*(lamda/600e-9)/lamda;
    beta = 2*pi*refractive_index/lamda;
    beta_0 = 2*pi/lamda;

    % 
    % figure(1)
    % imagesc(abs(ffted_kai(:,:,round(particle_center_index))));

    pixel_resolution = 2*pi/(2*k_max);
    k_resolution = k_max*2/(image_size-1);
    kai_resolution_k1 = 2*pi/(kai_size1*(particle_size1-1));
    kai_resolution_k2 = 2*pi/(kai_size1*(particle_size2-1));
    kai_resolution_k3 = 2*pi/(kai_size2*(particle_size3-1));
    Nor_cons1 = k_resolution/kai_resolution_k1;
    Nor_cons2 = k_resolution/kai_resolution_k2;
    Nor_cons3 = k_resolution/kai_resolution_k3;

    U_plus = zeros(image_size,image_size,phi_scan_number);
    U_minus = zeros(image_size,image_size,phi_scan_number);
    U_f = zeros(image_size,image_size,phi_scan_number);
    U_b = zeros(image_size,image_size,phi_scan_number);
    I_plus = zeros(image_size,image_size,phi_scan_number);
    I_minus = zeros(image_size,image_size,phi_scan_number);
    I_f = zeros(image_size,image_size,phi_scan_number);
    I_b = zeros(image_size,image_size,phi_scan_number);
    tic

% scanning phi-------------------------
    for m = 1:phi_scan_number
        tic
        phi = (m-1)*2*pi/phi_scan_number;
        Uplus_k = zeros(image_size,image_size,image_size);
        Nor_betaxin = beta*sin(theta)*sin(phi)/(kai_resolution_k1);
        Nor_betayin = beta*sin(theta)*cos(phi)/(kai_resolution_k2);
        Nor_gammain = beta*cos(theta)/(kai_resolution_k3);
        gammain = beta*cos(theta);
        betaxin = beta*sin(theta)*sin(phi);
        betayin = beta*sin(theta)*cos(phi);     




%getting U+ & U-    ----------------------------------------------

        for i = (1):(image_size)
            for j = (1):(image_size)
               for k = (1):(image_size)

                    if ((i-image_center_index)*k_resolution)^2+((j-image_center_index)*k_resolution)^2 +((k-image_center_index)*k_resolution)^2 >= (0.99*beta)^2
                        continue
                    end       
                    kx_index_f = round(Nor_cons1*(i-image_center_index)-Nor_betaxin) + particle_center_index1;
                    ky_index_f = round(Nor_cons2*(j-image_center_index)-Nor_betayin) + particle_center_index2; 
                    kz_index_f = round(Nor_cons3*(k-image_center_index)-Nor_gammain) + particle_center_index3;                      
                    k_perp_squ = ((i-image_center_index)*k_resolution)^2 +((j-image_center_index)*k_resolution)^2;
                    gamma = sqrt(beta^2 - k_perp_squ); 
                    kz_temp = (k-image_center_index)*k_resolution;            

                    Uplus_k(i,j,k) = (beta_0^2)*ffted_kai(kx_index_f,ky_index_f,kz_index_f)/((2*gamma)*(gamma-kz_temp));

               end

            end
        end
        Uplus_k = Uplus_k.*Aper;
%         for k = 1:image_size
%             Uplus_k(:,:,k) = Uplus_k(:,:,k).*Aper;
%         end
        
%     
%         U_space = ifftn(ifftshift(Uplus_k));
        U_space = fftshift(ifftn(ifftshift(Uplus_k)));
        U_plus(:,:,m) =U_space(:,:,image_center_index);
        I_plus(:,:,m) = abs(U_plus(:,:,m)).^2;
        toc
    end    

    for m = 1:phi_scan_number
        tic
        phi = (m-1)*2*pi/phi_scan_number;
        Uminus_k = zeros(image_size,image_size,image_size);
        Nor_betaxin = beta*sin(theta)*sin(phi)/(kai_resolution_k1);
        Nor_betayin = beta*sin(theta)*cos(phi)/(kai_resolution_k2);
        Nor_gammain = beta*cos(theta)/(kai_resolution_k3);
        gammain = beta*cos(theta);
        betaxin = beta*sin(theta)*sin(phi);
        betayin = beta*sin(theta)*cos(phi);     
 

        for i = (1):(image_size)
            for j = (1):(image_size)
                for k = (1):(image_size)

                    if ((i-image_center_index)*k_resolution)^2+((j-image_center_index)*k_resolution)^2+((k-image_center_index)*k_resolution)^2 >= (0.99*beta)^2
                        continue
                    end       
                    kx_index_b = round(Nor_cons1*(i-image_center_index)-Nor_betaxin) + particle_center_index1;
                    ky_index_b = round(Nor_cons2*(j-image_center_index)-Nor_betayin) + particle_center_index2; 
                    kz_index_b = round(Nor_cons3*(k-image_center_index)-Nor_gammain) + particle_center_index3;               
                    k_perp_squ = ((i-image_center_index)*k_resolution)^2 +((j-image_center_index)*k_resolution)^2;
                    gamma = sqrt(beta^2 - k_perp_squ); 
                    kz_temp = (k-image_center_index)*k_resolution;                 
                    Uminus_k(i,j,k) = (beta_0^2)*ffted_kai(kx_index_b,ky_index_b,kz_index_b)/((2*gamma)*(gamma+kz_temp));

                end
            end
        end
        Uminus_k = Uminus_k.*Aper;
%         for k = 1:image_size
%             Uminus_k(:,:,k) = Uminus_k(:,:,k).*Aper;
%         end
%     
%     
%         U_space = ifftn(ifftshift(Uminus_k));
        U_space = fftshift(ifftn(ifftshift(Uminus_k)));
        U_minus(:,:,m) =U_space(:,:,image_center_index);
        I_minus(:,:,m) = abs(U_minus(:,:,m)).^2;
        toc
    end


    for m = 1:phi_scan_number
        phi = (m-1)*2*pi/phi_scan_number;
        gammain = beta*cos(theta);
        betaxin = beta*sin(theta)*sin(phi);
        betayin = beta*sin(theta)*cos(phi); 
        for i = (1):(image_size)
            for j = (1):(image_size)
                  U_f(i,j,m) = U_plus(i,j,m) + alpha_f*exp(1i*gammain*z_plus*pixel_resolution+1i*betaxin*pixel_resolution*(i-image_center_index)+1i*betayin*pixel_resolution*(j-image_center_index));            
                  U_b(i,j,m) = U_minus(i,j,m) + alpha_b*exp(1i*gammain*z_minus*pixel_resolution+1i*betaxin*pixel_resolution*(i-image_center_index)+1i*betayin*pixel_resolution*(j-image_center_index));
            end
        end
        I_f(:,:,m) = abs(U_f(:,:,m)).^2;
        I_b(:,:,m) = abs(U_b(:,:,m)).^2;
    end

mp = mean(I_plus,3);
mm = mean(I_minus,3);
mf = mean(I_f,3);
mb =(mean(I_b,3));

fI_plus = fI_plus + spec(n,2)*mp*(1/lamda^2);
fI_minus = fI_minus + spec(n,2)*mm*(1/lamda^2);
fI_f = fI_f + spec(n,2)*mf*(1/lamda^2);
fI_b = fI_b + spec(n,2)*mb*(1/lamda^2);
fI_plus_sa = fI_plus_sa + spec(n,2)*I_plus(:,:,sa_index)*(1/lamda^2);
fI_minus_sa = fI_minus_sa + spec(n,2)*I_minus(:,:,sa_index)*(1/lamda^2);
fI_f_sa = fI_f_sa + spec(n,2)*I_f(:,:,sa_index)*(1/lamda^2);
fI_b_sa = fI_b_sa + spec(n,2)*I_b(:,:,sa_index)*(1/lamda^2);

%--------------------------------------------------------------------------
end


figure(1)
subplot(2,4,1);
imagesc(imresize(fI_plus_sa,0.5));
set(gca,'dataAspectRatio',[1,1,1])
set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
set(gca,'visible','off')
xlabel('DF forward scatter,phi=pi/6');
subplot(2,4,3);
imagesc(imresize(fI_minus_sa,0.5));
set(gca,'dataAspectRatio',[1,1,1])
set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
set(gca,'visible','off')
xlabel('DF back scatter,phi=pi/6');
subplot(2,4,5);
mp = mean(fI_plus,3);
imagesc(imresize(mp,0.5));
set(gca,'dataAspectRatio',[1,1,1])
set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
set(gca,'visible','off')
xlabel('DF forward scatter,scanning 120 phi');
subplot(2,4,7);
mm = mean(fI_minus,3);
imagesc(imresize(mm,0.5));
set(gca,'dataAspectRatio',[1,1,1])
set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
set(gca,'visible','off')
xlabel('DF back scatter,scanning 120 phi');
subplot(2,4,2);
imagesc(imresize(fI_f_sa,0.5));
set(gca,'dataAspectRatio',[1,1,1])
set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
set(gca,'visible','off')
xlabel('BF forward scatter,phi=pi/6');
subplot(2,4,4);
imagesc(imresize(fI_b_sa,0.5));
set(gca,'dataAspectRatio',[1,1,1])
set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
set(gca,'visible','off')
xlabel('BF back scatter,phi=pi/6');
subplot(2,4,6);
mf = mean(fI_f,3);
imagesc(imresize(mf,0.5));
set(gca,'dataAspectRatio',[1,1,1])
set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
set(gca,'visible','off')
xlabel('BF forward scatter,scanning 120 phi');
subplot(2,4,8);
mb =(mean(fI_b,3));
imagesc(imresize(mb,0.5));
set(gca,'dataAspectRatio',[1,1,1])
set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
set(gca,'visible','off')
xlabel('BF back scatter,scanning 120 phi');

figure(2)
for m = 1:phi_scan_number
    subplot(6,6,m);
    imagesc(imresize(I_plus(:,:,m),0.5));
    set(gca,'dataAspectRatio',[1,1,1])
    set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
    set(gca,'visible','off')
end


% for m = 1:phi_scan_number
%     I_f_new(:,:,m)=imresize(I_f(:,:,m),0.5);
% end

figure(3)
for m = 1:phi_scan_number
    subplot(6,6,m);
    imagesc(imresize(I_f(:,:,m),0.5));
    set(gca,'dataAspectRatio',[1,1,1])
    set(gca,'XTick',[], 'YTick', [], 'ZTick', [])
    set(gca,'visible','off')
end
