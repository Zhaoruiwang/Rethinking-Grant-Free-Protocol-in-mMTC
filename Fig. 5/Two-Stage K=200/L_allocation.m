function L_det = L_allocation(K_est)

if K_est>300
    L_det=round(1.1*K_est);
else
    if mod(K_est,10)==0
        value=K_est;
    else
        value=K_est-mod(K_est,10)+10;
    end
    switch value
        case 10
            L_det=15;
        case 20
            L_det=25;
        case 30
            L_det=35;
        case 40
            L_det=45;
        case 50
            L_det=55;
        case 60
            L_det=65;
        case 70
            L_det=75;
        case 80
            L_det=76;
        case 90
            L_det=77;
        case 100
            L_det=78;
        case 110
            L_det=80;
        case 120
            L_det=90;
        case 130
            L_det=100;
        case 140
            L_det=110;
        case 150
            L_det=120;
        case 160
            L_det=130;
        case 170
            L_det=140;
        case 180
            L_det=150;
        case 190
            L_det=160;
        case 200
            L_det=170;
        case 210
            L_det=180;
        case 220
            L_det=190;
        case 230
            L_det=200;
        case 240
            L_det=210;
        case 250
            L_det=220;
        case 260
            L_det=230;
        case 270
            L_det=240;
        case 280
            L_det=250;
        case 290
            L_det=260;
        case 300
            L_det=270;
    end
end
end