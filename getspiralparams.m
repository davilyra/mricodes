function [kfile,spatfov,spatres,nintl,nread] = getspiralparams(spiralid)

nintl=[];
nread=[];

switch(spiralid)
    case 0,
        kfile='recon65mm20cm'; spatfov = 200; spatres = 6.5;
    case 2,
        kfile='recon25mm20cm'; spatfov = 200; spatres = 2.5;
    case 5,
        kfile='recon7mm25cm'; spatfov = 250; spatres = 7;
    case 6,
        kfile='recon2mm15cm'; spatfov = 150; spatres = 2;
    case 7,
        kfile='recon7mm25cm4vd'; spatfov = 250; spatres = 7;
    case 8,
        kfile='recon47mm25cm4vd'; spatfov = 250; spatres = 4.7;
    case 9,
        kfile='recon47mm25cm4vd2intl'; spatfov = 250; spatres = 4.7;
    case 10,
        kfile='recon36mm25cm4vd3intl'; spatfov = 250; spatres = 3.6;
    case 11,
        kfile='recon3mm25cm4vd4intl'; spatfov = 250; spatres = 3; 
    case 12,
        kfile='recon36mm25cm1vd6intl'; spatfov = 250; spatres = 3.6; 
    case 13,
        kfile='recon47mm125cm4vd1intl'; spatfov = 125; spatres = 4.7;         
    case 14,
        kfile='recon225mm25cm4vd6intl'; spatfov = 250; spatres = 2.25; nintl = 6; nread = 1028;
    case 15,
        kfile='recon225mm25cm1vd12intl'; spatfov = 250; spatres = 2.25; nintl = 12; nread = 972;       
    case 16,
        kfile='recon16cm09mm16intl4vd'; spatfov = 160; spatres = 0.9;                
    case 17,
        kfile='recon10cm068mm16intl4vd'; spatfov = 100; spatres = 0.68;           
    case 18,
        kfile='recon8cm1mm7intl4vd'; spatfov = 80; spatres = 1; nintl = 7; nread = 970;
    case 19,
        kfile='recon16cm14mm8intl4vd'; spatfov = 160; spatres = 1.4; nintl = 8; nread = 1012;        
    case 20,
        kfile='recon13cm3mm2intl4vd'; spatfov = 130; spatres = 3;
    case 21,
        kfile='recon14cm2mm4intl4vd'; spatfov = 140; spatres = 2;
    case 22,
        kfile='recon35cm14mm12intl6vd'; spatfov = 350; spatres = 1.4;        
    case 23,
        kfile='recon35cm065mm48intl6vd'; spatfov = 350; spatres = 0.65;                
    case 24,
        kfile='recon13cm033mm150intl1vd'; spatfov = 130; spatres = 0.33;                                
    case 25,
        kfile='recon16cm1mm24intl1vd'; spatfov = 160; spatres = 1;  nintl = 24; nread = 1048;                                       
    otherwise,
        error(['unexpected spiralid value: ',num2str(spiralid)]);
end;