close all
clc
clear

myPSNR=[];
mySSIM=[];
for gausigma=[10 20 30 40 50]
    
        picindex = {'1B.bmp','1J.jpg','1P.png','2P.png','3P.png','4P.png','5P.png','6P.png','7P.png','8P.png','9P.png','10P.png','1T.tif','2T.tif','1TF.tiff','2TF.tiff','3TF.tiff','4TF.tiff','5TF.tiff','6TF.tiff'};

    for i=1:length(picindex)
        
        fprintf('the %d th picture\n',i);
        str=['..\Images\',picindex{i}];
        ima=imread(str);
        if length(size(ima))>2
            ima=rgb2gray(ima);
        end
        
        ima=imresize(ima,[256 256]);
        [m,n]=size(ima);

        str=['..\results\',picindex{i},'_.tiff'];
        imwrite(uint8(ima),str);
        ima=double(ima);
        Original_I=ima;
        
        
        Input=ima+gausigma*randn(size(ima));%mixed noise
        str=['..\results\',picindex{i},'_GS',num2str(gausigma),'_.tiff'];
        imwrite(uint8(Input),str);
        
        InputI = Input;
        
        
        %MES
        MES_Input = sum(sum((double(InputI)-double(Original_I)).*(double(InputI)-double(Original_I))))/size(Original_I,1)/size(Original_I,2);
        %PSNR
        PSNR_Input = 10 * log10(255 * 255/ MES_Input)
        
        SSIM_Input = ssim(ima,InputI)
        
        
        
        
        for  Cluster=50
            
            for SS=50
                
                pcadim=6;
                S=8;
                K=3;
                
                se=strel('square',3);
                I=InputI;
                
                sigma=gausigma;
                
                W0=I;
                
                [MAP]=LNLMMMAP(W0,sigma,Cluster,SS,S,K,pcadim);
                
                %LNLMM CO-Closing
                I_DC=imdilate(I,se);
                I_DN=NLMGNNDilation(I,MAP);
                I_LND=max(I_DC,I_DN);
                
                I_CC=imerode(I_LND,se);
                I_CN=NLMGNNErosion(I_LND,MAP);
                I_LNC=min(I_CC,I_CN);
                
                
                %LNLMM CO
                I_CEC=imerode(I_LNC,se);
                I_CEN=NLMGNNErosion(I_LNC,MAP);
                I_LNCE=min(I_CEC,I_CEN);
                
                I_COC=imdilate(I_LNCE,se);
                I_CON=NLMGNNDilation(I_LNCE,MAP);
                R_LNCO=max(I_COC,I_CON);
                
                %LNLMM OC-Opening
                I_EC=imerode(I,se);
                I_EN=NLMGNNErosion(I,MAP);
                I_LNE=min(I_EC,I_EN);
                
                I_OC=imdilate(I_LNE,se);
                I_ON=NLMGNNDilation(I_LNE,MAP);
                I_LNO=max(I_OC,I_ON);
                
                %LNLMM OC
                I_ODC=imdilate(I_LNO,se);
                I_ODN=NLMGNNDilation(I_LNO,MAP);
                I_LNOD=max(I_ODC,I_ODN);
                
                I_OCC=imerode(I_LNOD,se);
                I_OCN=NLMGNNErosion(I_LNOD,MAP);
                R_LNOC=min(I_OCC,I_OCN);
                
                
                LNLMMReconstruction=(double(R_LNCO)+double(R_LNOC))/2;
                
                figure,imshow(uint8(LNLMMReconstruction));title('LNLMMReconstruction');
                str=['..\results\',picindex{i},'_LNLMM',num2str(SS),num2str(gausigma),'_.tiff'];
                imwrite(uint8(LNLMMReconstruction),str);
                
                
                MES_LNLMMReconstruction= sum(sum((double(LNLMMReconstruction)-double(Original_I)).*(double(LNLMMReconstruction)-double(Original_I))))/size(Original_I,1)/size(Original_I,2);
                
                PSNR_LNLMMReconstruction = 10 * log10(255 * 255/ MES_LNLMMReconstruction)
                SSIM_LNLMMReconstruction = ssim(ima,LNLMMReconstruction)
                
            end
            
        end
        
        myPSNR=[myPSNR; i gausigma  PSNR_Input  PSNR_LNLMMReconstruction];
        mySSIM=[mySSIM; i gausigma  SSIM_Input SSIM_LNLMMReconstruction];
        
    end
    
    close all
    
    
end
save('MyPSNRData', 'myPSNR');
save('MySSIMData', 'mySSIM');



