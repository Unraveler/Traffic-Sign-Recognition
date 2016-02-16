%% Programa protótipo para deteção e georeferenciação de sinalização vertical
% -------------------------------------------------------------------------
clc;
clear all;
close all;
format long g;
% Projeto em Informação Geográfica 2012/2013
% Mestrado em Engenharia Geográfica
% Autor: Rui Jorge Abrunhosa Nunes nº32092
% -------------------------------------------------------------------------
%% 1ª Parte: Inputs
    % Abertura do ficheiro para escrita
    centroides = fopen('centroides.txt','w');
    resultados = fopen('resultados.txt','w');
    
    %Leitura das colunas do ficheiro de orientações externas. Guardar em
    %vector
    load ori_externas.txt;
    X0 = ori_externas(:,2);
    Y0 = ori_externas(:,3);
    Z0 = ori_externas(:,4);
    heading = ori_externas(:,5);
    pitch = ori_externas(:,6);
    roll = ori_externas(:,7);
    sigma_heading = ori_externas(:,8);
    sigma_pitch = ori_externas(:,9);
    sigma_roll = ori_externas(:,10);
    frame_num = ori_externas(:,1);
    
    %Tamnho de um pixel
    pixX = 0.00465;
    pixY = pixX;
   
    % Abertura do video
    video = mmreader('video.avi'); %Leitura do ficheiro de video
    num_frames = video.NumberOfFrames; %Variavel que contem o numero de frames do video

    %Ciclo para tratar independentemene cada frame
    for i = 1:num_frames
        % Leitura das frames
        frame = read(video,i);
        realframe = i+579;
        %% 2ª Parte: Processamento Digital de Imagem
        % Segmentação de cor
            %(Conversão RGB para HSV)
            frameHSV = rgb2hsv(frame);
            %Compoenntes da imagem
            H = frameHSV(:,:,1);
            S = frameHSV(:,:,2);
            V = frameHSV(:,:,3);
        %Threshold de cada componente
        Thresh_H = im2bw(H,0.5);
        Thresh_S = im2bw(S,0.5);
        Thresh_V = im2bw(V);
        %Filtrar as componentes em Threshold
            %Preencher áreas abertas
            Thresh_H = imfill(Thresh_H,'holes');
            Thresh_S = imfill(Thresh_S,'holes');
            Thresh_V = imfill(Thresh_V,'holes');
            %Multiplicar componente H por S
            Thresh = Thresh_H .* Thresh_S;
            %Eliminar pixeis que não se encontrem conectados
            Thresh = bwareaopen(Thresh,1000,8);

            figure(1),imshow(frame), hold on;
            
            % Apagar todos os objectos que não são sinais da imagem
            frame_sinais = frame; %Criação de uma matriz auxiliar
            for m = 1: size(frame,1)
                for n = 1: size(frame,2)
                    if Thresh(m,n) == 0
                        frame_sinais(m,n,:) = 256; %O que não é sinal, fica a branco
                    else
                        frame_sinais(m,n,:) = frame(m,n,:); %O que é sinal, mantém a aparência original
                    end
                end
            end
         
        %Detetar sinais e correlacionar
        ObjectosConectados = bwlabel(Thresh); %Esta função serve para marcar todos os objectos que estão presentes na imagem binária
        Propriedades = regionprops(ObjectosConectados,'Centroid','Eccentricity','FilledArea','Perimeter','BoundingBox'); %Esta função serve para detetar certas caracteristicas dos objectos
           %Ciclo para marcar os sinais detetados (se não for feito, só é marcado um sinal, caso existam 2 ou mais)
            for x = 1:size(Propriedades,1)
                if((Propriedades(x).Eccentricity >=0.05) && (Propriedades(x).Eccentricity <=0.45))
                    AreaCheia = Propriedades(x).FilledArea(1); %Calculo da área do sinal
                    Raio = round(sqrt(AreaCheia/pi())); %Calculo do raio da figura para o corte
                    Raio = Raio + 3;
                    CentroX = round(Propriedades(x).Centroid(1)); %Guardar a coordenada X do centroide do sinal
                    CentroY = round(Propriedades(x).Centroid(2)); %Guardar a coordenada Y do centroide do sinal
                    plot(CentroX,CentroY,'r*');
                end
            end
            
            % Ciclo para destacar os sinais e fazer correlação
            for x = 1:size(Propriedades,1)
                if((Propriedades(x).Eccentricity(1) >=0.05) && (Propriedades(x).Eccentricity(1) <=0.45)) %Filtragem pela forma do objecto. Se o objecto apresentar valores de forma entre estes limites, então é um sinal redondo ou triangular
                    CentroX = round(Propriedades(x).Centroid(1)); %Guardar a coordenada X do centroide do sinal
                    CentroY = round(Propriedades(x).Centroid(2)); %Guardar a coordenada Y do centroide do sinal
                    AreaCheia = Propriedades(x).FilledArea(1); %Calculo da área do sinal
                    Raio = round(sqrt(AreaCheia/pi())); %Calculo do raio da figura para o corte
                    Raio = Raio + 4;
                   
                    %Definir canto superior esquerdo e inferior direito,
                    %para destacar o sinal.
                    Canto1X = CentroX-Raio;
                    Canto1Y = CentroY-Raio;
                    Canto4X = CentroX+Raio;
                    Canto4Y = CentroY+Raio;
                    %Define uma margem para o corte do sinal
                    margem = 0;
                    %Corta a volta do sinal e cria uma nova matriz
                    corte = imcrop(frame,[Canto1X-margem Canto1Y-margem (Canto4X-Canto1X)+3*margem (Canto4Y-Canto1Y)+3*margem]); % Fazer um corte na imagem de modo a criar uma pequena matriz só com o sinal
                    corte_sinais = imcrop(frame_sinais,[Canto1X-margem Canto1Y-margem (Canto4X-Canto1X)+3*margem (Canto4Y-Canto1Y)+3*margem]); % Fazer um corte na imagem de modo a criar uma pequena matriz só com o sinal
                    %Escalar os sinais cortados
                    corte =  imresize(corte,[100,100],'cubic');
                    corte_sinais = imresize(corte_sinais,[100,100],'cubic');
                    %Correlação com modelos
                        % Abertura dos modelos
                        template00 = imread('PerigoCurvaDireita.bmp'); %Sinal de curva perigosa a direita
                        template01 = imread('Entroncamento.bmp'); %Sinal de perigo via com entroncamento à direita
                        template02 = imread('Prioridade.bmp'); %Sinal de cedência de prioridade
                        template03 = imread('ObrigatorioFrente.bmp'); %Sinal de obrigatório seguir em fente
                        template04 = imread('ProibidoEsquerda.bmp'); %Sinal de proibido virar à esquerda
                        template05 = imread('ProibidoUltrapassar.bmp'); %Sinal de proibição de ultrapassagem
                        template06 = imread('Proibido40.bmp'); %Sinal de proibição de velocidade > 40km
                        template07 = imread('Proibido60.bmp'); %Sinal de proibição de velocidade > 60km
                        template08 = imread('Proibido80.bmp'); %Sinal de proibição de velocidade > 80km
            
                        % Escalar os modelos
                        template00 = imresize(template00,[100,100],'cubic');
                        template01 = imresize(template01,[100,100],'cubic');
                        template02 = imresize(template02,[100,100],'cubic');
                        template03 = imresize(template03,[100,100],'cubic');
                        template04 = imresize(template04,[100,100],'cubic');
                        template05 = imresize(template05,[100,100],'cubic');
                        template06 = imresize(template06,[100,100],'cubic');
                        template07 = imresize(template07,[100,100],'cubic');
                        template08 = imresize(template08,[100,100],'cubic');
                     
                        % Correlação de imagens
                        corr00 = normxcorr2(template00(:,:,2),corte_sinais(:,:,2));
                        corr01 = normxcorr2(template01(:,:,2),corte_sinais(:,:,2));
                        corr02 = normxcorr2(template02(:,:,2),corte_sinais(:,:,2));
                        corr03 = normxcorr2(template03(:,:,1),corte_sinais(:,:,1));
                        corr04 = normxcorr2(template04(:,:,2),corte_sinais(:,:,2));
                        corr05 = normxcorr2(template05(:,:,2),corte_sinais(:,:,2));
                        corr06 = normxcorr2(template06(:,:,2),corte_sinais(:,:,2));
                        corr07 = normxcorr2(template07(:,:,2),corte_sinais(:,:,2));
                        corr08 = normxcorr2(template08(:,:,3),corte_sinais(:,:,3));

                        max_corr00 = max(max(corr00(:)));
                        max_corr01 = max(max(corr01(:)));
                        max_corr02 = max(max(corr02(:)));
                        max_corr03 = max(max(corr03(:)));
                        max_corr04 = max(max(corr04(:)));
                        max_corr05 = max(max(corr05(:)));
                        max_corr06 = max(max(corr06(:)));
                        max_corr07 = max(max(corr07(:)));
                        max_corr08 = max(max(corr08(:)));
                        
                         %Coordenadas pixel referidas ao centro da imagem
                            CentroX_CF = CentroY - ((size(frame,2))/2);
                            CentroY_CF = ((size(frame,1))/2) - CentroX;
                         %Coordenadas foto para qualquer ponto
                            CentroX = CentroX_CF*pixX;
                            CentroY = CentroY_CF*pixY;
                    
                        if(max_corr00>=0.3164 && max_corr00 <=0.3166 || max_corr00>=0.4151 && max_corr00 <=0.4153)
                           matriz(i,1) = 100;
                           matriz(i,2) = realframe;
                           matriz(i,3) = ori_externas(i,1);
                           matriz(i,4) = CentroX;
                           matriz(i,5) = CentroY;
                           matriz(i,6) = ori_externas(i,2);
                           matriz(i,7) = ori_externas(i,3);
                           matriz(i,8) = ori_externas(i,4);
                           matriz(i,9) = ori_externas(i,5);
                           matriz(i,10) = ori_externas(i,6);
                           matriz(i,11) = ori_externas(i,7);
                           matriz(i,12) = ori_externas(i,8);
                           matriz(i,13) = ori_externas(i,9);
                           matriz(i,14) = ori_externas(i,10);
                        elseif(max_corr01>=0.2857 && max_corr01 <=0.2859 || max_corr01>=0.2897 && max_corr01 <=0.2899)
                            matriz(i,1) = 101;
                            matriz(i,2) = realframe;
                            matriz(i,3) = ori_externas(i,1);
                            matriz(i,4) = CentroX;
                            matriz(i,5) = CentroY;
                            matriz(i,6) = ori_externas(i,2);
                            matriz(i,7) = ori_externas(i,3);
                            matriz(i,8) = ori_externas(i,4);
                            matriz(i,9) = ori_externas(i,5);
                            matriz(i,10) = ori_externas(i,6);
                            matriz(i,11) = ori_externas(i,7);
                            matriz(i,12) = ori_externas(i,8);
                            matriz(i,13) = ori_externas(i,9);
                            matriz(i,14) = ori_externas(i,10);  
                        elseif(max_corr02>= 0.2797 && max_corr02 <=0.2799 || max_corr02>=0.2645 && max_corr02 <=0.2647)
                            matriz(i,1) = 102;
                            matriz(i,2) = realframe;
                            matriz(i,3) = ori_externas(i,1);
                            matriz(i,4) = CentroX;
                            matriz(i,5) = CentroY;
                            matriz(i,6) = ori_externas(i,2);
                            matriz(i,7) = ori_externas(i,3);
                            matriz(i,8) = ori_externas(i,4);
                            matriz(i,9) = ori_externas(i,5);
                            matriz(i,10) = ori_externas(i,6);
                            matriz(i,11) = ori_externas(i,7);
                            matriz(i,12) = ori_externas(i,8);
                            matriz(i,13) = ori_externas(i,9);
                            matriz(i,14) = ori_externas(i,10);
                        elseif(max_corr03>=0.7646 && max_corr03 <=0.7648 || max_corr03>=0.8793 && max_corr03 <=0.8795)
                            matriz(i,1) = 103;
                            matriz(i,2) = realframe;
                            matriz(i,3) = ori_externas(i,1);
                            matriz(i,4) = CentroX;
                            matriz(i,5) = CentroY;
                            matriz(i,6) = ori_externas(i,2);
                            matriz(i,7) = ori_externas(i,3);
                            matriz(i,8) = ori_externas(i,4);
                            matriz(i,9) = ori_externas(i,5);
                            matriz(i,10) = ori_externas(i,6);
                            matriz(i,11) = ori_externas(i,7);
                            matriz(i,12) = ori_externas(i,8);
                            matriz(i,13) = ori_externas(i,9);
                            matriz(i,14) = ori_externas(i,10);
                        elseif(max_corr04>=0.5697 && max_corr04 <=0.5699 || max_corr04>=0.7957 && max_corr04 <=0.7959)
                            matriz(i,1) = 104;
                            matriz(i,2) = realframe;
                            matriz(i,3) = ori_externas(i,1);
                            matriz(i,4) = CentroX;
                            matriz(i,5) = CentroY;
                            matriz(i,6) = ori_externas(i,2);
                            matriz(i,7) = ori_externas(i,3);
                            matriz(i,8) = ori_externas(i,4);
                            matriz(i,9) = ori_externas(i,5);
                            matriz(i,10) = ori_externas(i,6);
                            matriz(i,11) = ori_externas(i,7);
                            matriz(i,12) = ori_externas(i,8);
                            matriz(i,13) = ori_externas(i,9);
                            matriz(i,14) = ori_externas(i,10);
                        elseif(max_corr05>=0.5016 && max_corr05 <=0.5018 || max_corr05>=0.6459 && max_corr05 <=0.6461)
                            matriz(i,1) = 105;
                            matriz(i,2) = realframe;
                            matriz(i,3) = ori_externas(i,1);
                            matriz(i,4) = CentroX;
                            matriz(i,5) = CentroY;
                            matriz(i,6) = ori_externas(i,2);
                            matriz(i,7) = ori_externas(i,3);
                            matriz(i,8) = ori_externas(i,4);
                            matriz(i,9) = ori_externas(i,5);
                            matriz(i,10) = ori_externas(i,6);
                            matriz(i,11) = ori_externas(i,7);
                            matriz(i,12) = ori_externas(i,8);
                            matriz(i,13) = ori_externas(i,9);
                            matriz(i,14) = ori_externas(i,10);
                        elseif(max_corr06>=0.5281 && max_corr06 <=0.5283 || max_corr06>=0.7064 && max_corr06 <=0.7066)
                            matriz(i,1) = 106;
                            matriz(i,2) = realframe;
                            matriz(i,3) = ori_externas(i,1);
                            matriz(i,4) = CentroX;
                            matriz(i,5) = CentroY;
                            matriz(i,6) = ori_externas(i,2);
                            matriz(i,7) = ori_externas(i,3);
                            matriz(i,8) = ori_externas(i,4);
                            matriz(i,9) = ori_externas(i,5);
                            matriz(i,10) = ori_externas(i,6);
                            matriz(i,11) = ori_externas(i,7);
                            matriz(i,12) = ori_externas(i,8);
                            matriz(i,13) = ori_externas(i,9);
                            matriz(i,14) = ori_externas(i,10);
                        elseif(max_corr07>=0.5186 && max_corr07 <=0.5188 || max_corr07>=0.5676 && max_corr07 <=0.5678)
                            matriz(i,1) = 107;
                            matriz(i,2) = realframe;
                            matriz(i,3) = ori_externas(i,1);
                            matriz(i,4) = CentroX;
                            matriz(i,5) = CentroY;
                            matriz(i,6) = ori_externas(i,2);
                            matriz(i,7) = ori_externas(i,3);
                            matriz(i,8) = ori_externas(i,4);
                            matriz(i,9) = ori_externas(i,5);
                            matriz(i,10) = ori_externas(i,6);
                            matriz(i,11) = ori_externas(i,7);
                            matriz(i,12) = ori_externas(i,8);
                            matriz(i,13) = ori_externas(i,9);
                            matriz(i,14) = ori_externas(i,10);
                        elseif(max_corr08>=0.5879 && max_corr08 <=0.5881 || max_corr08>=0.6304 && max_corr08 <=0.6306)
                            matriz(i,1) = 108;
                            matriz(i,2) = realframe;
                            matriz(i,3) = ori_externas(i,1);
                            matriz(i,4) = CentroX;
                            matriz(i,5) = CentroY;
                            matriz(i,6) = ori_externas(i,2);
                            matriz(i,7) = ori_externas(i,3);
                            matriz(i,8) = ori_externas(i,4);
                            matriz(i,9) = ori_externas(i,5);
                            matriz(i,10) = ori_externas(i,6);
                            matriz(i,11) = ori_externas(i,7);
                            matriz(i,12) = ori_externas(i,8);
                            matriz(i,13) = ori_externas(i,9);
                            matriz(i,14) = ori_externas(i,10);
                        end
               end
            end
           
    end
matriz = sortrows(matriz,1); %Ordena as colunas da matriz por identificação do sinal de transito

%Ciclo para criar uma matriz temporária com todos os zeros que aparecem
%indevidamente antes do que interessa
for m = 1:size(matriz,1)
    for n = 1:size(matriz,2)
        if(matriz(m,n) == 0)
            matriz_temp(m,n) = matriz(m,n);
        end
    end
end

matriz2 = matriz(size(matriz_temp,1)+1:size(matriz,1),1:size(matriz,2)); %Corta os zeros da matriz anterior e cria uma apenas com os valores de interesse
dlmwrite('centroides.txt', matriz2,'delimiter','\t','precision','%7.6f'); %Escreve a matriz em ficheiro. Colunas delimitadas por tabulação
fclose(centroides);
    %% 3ª Parte: Georreferenciação
    
    load centroides.txt;

%Criar vectores independentes para cada variável

    %Identificações
    id_sinal = centroides(:,1);
    %Coordenadas dos pontos nas fotos
    centrox = centroides(:,4);
    centroy = centroides(:,5);
    %Parâmetros externos
    X0 = centroides(:,6);
    Y0 = centroides(:,7);
    Z0 = centroides(:,8);
    heading = centroides(:,9);
    pitch = centroides(:,10);
    roll = centroides(:,11);
    %Parâmetros internos
    x0 = 0.001086*10e-3;
    y0 = 0.005296*10e-3;
    c = 6.089964*10e-3;
    
    %Separar por fotos
    %Foto1
    for m = 1: size(centrox,1)/2
        id_sinal1(m,1) = id_sinal(2*m-1,1);
        x1(m,1) = centrox(2*m-1,1);
        y1(m,1) = centroy(2*m-1,1);
        X0_foto1(m,1) = X0(2*m-1,1);
        Y0_foto1(m,1) = Y0(2*m-1,1);
        Z0_foto1(m,1) = Z0(2*m-1,1);
        heading_foto1(m,1) = heading(2*m-1,1);
        pitch_foto1(m,1) = pitch(2*m-1,1);
        roll_foto1(m,1) = roll(2*m-1,1);
    end   
    
    %Foto2
    for m = 1: size(centrox,1)/2
        id_sinal2(m,1) = id_sinal(2*m,1);
        x2(m,1) = centrox(2*m,1);
        y2(m,1) = centroy(2*m,1);
        X0_foto2(m,1) = X0(2*m,1);
        Y0_foto2(m,1) = Y0(2*m,1);
        Z0_foto2(m,1) = Z0(2*m,1);
        heading_foto2(m,1) = heading(2*m,1);
        pitch_foto2(m,1) = pitch(2*m,1);
        roll_foto2(m,1) = roll(2*m,1);
    end
%--------------------------------------------------------------------------

 %Escrita do cabeçalho do ficheiro de resultados
    fprintf(resultados,'Ponto \t\t\t\t\t\t\t\t\tM\t\t\tP\t\t\tH\n');
    fprintf(resultados,'-----------------------------------------------------------------------\n');
    
    for k = 1:size(id_sinal1,1)
%% Calculo das Coordenadas Terreno segundo as Equações de Colineariedade---
    %Conversor de grados para garus
    conv = (pi/180);
    %Matriz de Rotação para a Fotografia 1
    Rot_foto1(1,1) = cos(pitch_foto1(k)*conv)*cos(heading_foto1(k)*conv);
    Rot_foto1(1,2) = cos(roll_foto1(k)*conv)*sin(heading_foto1(k)*conv);
    Rot_foto1(1,3) = sin(pitch_foto1(k)*conv)*cos(heading_foto1(k)*conv) - sin(heading_foto1(k)*conv)*sin(roll_foto1(k)*conv);
    Rot_foto1(2,1) = -cos(pitch_foto1(k)*conv)*sin(heading_foto1(k)*conv);
    Rot_foto1(2,2) = cos(roll_foto1(k)*conv)*cos(heading_foto1(k)*conv);
    Rot_foto1(2,3) = sin(heading_foto1(k)*conv)*sin(pitch_foto1(k)*conv) + sin(roll_foto1(k)*conv)*cos(heading_foto1(k)*conv);
    Rot_foto1(3,1) = sin(pitch_foto1(k)*conv);
    Rot_foto1(3,2) = sin(roll_foto1(k)*conv);
    Rot_foto1(3,3) = cos(roll_foto1(k)*conv)*cos(pitch_foto1(k)*conv);
    
    %Matriz de Rotação A Fotografia 2
    Rot_foto2(1,1) = cos(pitch_foto2(k)*conv)*cos(heading_foto2(k)*conv);
    Rot_foto2(1,2) = cos(roll_foto2(k)*conv)*sin(heading_foto2(k)*conv);
    Rot_foto2(1,3) = sin(pitch_foto2(k)*conv)*cos(heading_foto2(k)*conv) - sin(heading_foto2(k)*conv)*sin(roll_foto2(k)*conv);
    Rot_foto2(2,1) = -cos(pitch_foto2(k)*conv)*sin(heading_foto2(k)*conv);
    Rot_foto2(2,2) = cos(roll_foto2(k)*conv)*cos(heading_foto2(k)*conv);
    Rot_foto2(2,3) = sin(heading_foto2(k)*conv)*sin(pitch_foto2(k)*conv) + sin(roll_foto2(k)*conv)*cos(heading_foto2(k)*conv);
    Rot_foto2(3,1) = sin(pitch_foto2(k)*conv);
    Rot_foto2(3,2) = sin(roll_foto2(k)*conv);
    Rot_foto2(3,3) = cos(roll_foto2(k)*conv)*cos(pitch_foto2(k)*conv);

        %CÁLCULO DAS APROXIMAÇÕES INICIAIS
        
        %Simplificação das equações
        %Fotografia 1
        Kx_foto1 = (Rot_foto1(1,1)*(x1(k)-x0) + Rot_foto1(1,2)*(y1(k)-y0) - Rot_foto1(1,3)*c) / (Rot_foto1(3,1)*(x1(k)-x0) + Rot_foto1(3,2)*(y1(k)-y0) - Rot_foto1(3,3)*c); 
        Ky_foto1 = (Rot_foto1(2,1)*(x1(k)-x0) + Rot_foto1(2,2)*(y1(k)-y0) - Rot_foto1(2,3)*c) / (Rot_foto1(3,1)*(x1(k)-x0) + Rot_foto1(3,2)*(y1(k)-y0) - Rot_foto1(3,3)*c);
        %Fotografia 2
        Kx_foto2 = (Rot_foto2(1,1)*(x2(k)-x0) + Rot_foto2(1,2)*(y2(k)-y0) - Rot_foto2(1,3)*c) / (Rot_foto2(3,1)*(x2(k)-x0) + Rot_foto2(3,2)*(y2(k)-y0) - Rot_foto2(3,3)*c);
        Ky_foto2 = (Rot_foto2(3,1)*(x2(k)-x0) + Rot_foto2(3,2)*(y2(k)-y0) - Rot_foto2(3,3)*c) / (Rot_foto2(3,1)*(x2(k)-x0) + Rot_foto2(3,2)*(y2(k)-y0) - Rot_foto2(3,3)*c);
     
        Z = 0;
        
        %Cálculo das coordenadas dos Pontos no terreno, da Fotografia 1
        X_foto1(k) = X0_foto1(k) + (Z-Z0_foto1(k)) * Kx_foto1;
        Y_foto1(k) = Y0_foto1(k) + (Z-Z0_foto1(k)) * Ky_foto1;
        
        %Cálculo das coordenadas dos Pontos no terreno, da Fotografia 1
        X_foto2(k) = X0_foto2(k) + (Z-Z0_foto2(k)) * Kx_foto2;
        Y_foto2(k) = Y0_foto2(k) + (Z-Z0_foto2(k)) * Ky_foto2;
        
       
%--------------------------------------------------------------------------
%% Ajustamento Paramétrico Não Linear
%Número total de observações
n = 4;

%Número de observações necessárias
n0 = 3;

%Graus de Liberdade
df = n-n0;

%Variância à priori
sigma0 = 1;

%Matriz de Variâncias/Covariâncias das observações
Cl = eye(n);

%Matriz dos pesos da observações
Pl = sigma0*inv(Cl);

%Vector de Observações
l = [x1(k);y1(k);x2(k);y2(k)];        

%Vector de parâmetros desconhecidos (com aproximações inicias)
X_desc = (X_foto1(k) + X_foto2(k))/2;
Y_desc = (Y_foto1(k) + Y_foto2(k))/2;
Z_desc = (X0_foto2(k)-Z0_foto2(k) * Kx_foto2+Z0_foto1(k)*Kx_foto1-X0_foto1(k)) / (Kx_foto1 - Kx_foto2);

resultado = [X_desc; Y_desc; Z_desc];

%----------------------------INICIO DO 1ºCLICLO------------------------------
precisao = 1e-4;
delta = 1;
interaccao = 0;
while max(abs(delta)) > precisao
    
    %SIMPLIFICAÇÃO DO MODELO MATEMÁTICO
    %Fotografia 1
    Nx_foto1 = Rot_foto1(1,1)*(resultado(1)-X0_foto1(k)) + Rot_foto1(2,1)*(resultado(2)-Y0_foto1(k)) + Rot_foto1(3,1)*(resultado(3)-Z0_foto1(k));
    Ny_foto1 = Rot_foto1(1,2)*(resultado(1)-X0_foto1(k)) + Rot_foto1(2,2)*(resultado(2)-Y0_foto1(k)) + Rot_foto1(3,2)*(resultado(3)-Z0_foto1(k));
    D_foto1 =  Rot_foto1(1,3)*(resultado(1)-X0_foto1(k)) + Rot_foto1(2,3)*(resultado(2)-Y0_foto1(k)) + Rot_foto1(3,3)*(resultado(3)-Z0_foto1(k));
    
    %Fotografia 2
    Nx_foto2 = Rot_foto2(1,1)*(resultado(1)-X0_foto2(k)) + Rot_foto2(2,1)*(resultado(2)-Y0_foto2(k)) + Rot_foto2(3,1)*(resultado(3)-Z0_foto2(k));
    Ny_foto2 = Rot_foto2(1,2)*(resultado(1)-X0_foto2(k)) + Rot_foto2(2,2)*(resultado(2)-Y0_foto2(k)) + Rot_foto2(3,2)*(resultado(3)-Z0_foto2(k));
    D_foto2 =  Rot_foto2(1,3)*(resultado(1)-X0_foto2(k)) + Rot_foto2(2,3)*(resultado(2)-Y0_foto2(k)) + Rot_foto2(3,3)*(resultado(3)-Z0_foto2(k));
    
    %Matriz dos coeficientes das correcções aos parâmetros desconhecidos (Primeira
    %Matriz de Configuração)
    A(1,1) = (c/D_foto1^2)*(-D_foto1*Rot_foto1(1,1) + Nx_foto1*Rot_foto1(1,3));
    A(1,2) = (c/D_foto1^2)*(-D_foto1*Rot_foto1(2,1) + Nx_foto1*Rot_foto1(2,3));
    A(1,3) = (c/D_foto1^2)*(-D_foto1*Rot_foto1(3,1) + Nx_foto1*Rot_foto1(3,3));
    A(2,1) = (c/D_foto1^2)*(-D_foto1*Rot_foto1(1,2) + Ny_foto1*Rot_foto1(1,3));
    A(2,2) = (c/D_foto1^2)*(-D_foto1*Rot_foto1(2,2) + Ny_foto1*Rot_foto1(2,3));
    A(2,3) = (c/D_foto1^2)*(-D_foto1*Rot_foto1(3,2) + Ny_foto1*Rot_foto1(3,3));
    A(3,1) = (c/D_foto2^2)*(-D_foto2*Rot_foto2(1,1) + Nx_foto2*Rot_foto2(1,3));
    A(3,2) = (c/D_foto2^2)*(-D_foto2*Rot_foto2(2,1) + Nx_foto2*Rot_foto2(2,3));
    A(3,3) = (c/D_foto2^2)*(-D_foto2*Rot_foto2(3,1) + Nx_foto2*Rot_foto2(3,3));
    A(4,1) = (c/D_foto2^2)*(-D_foto2*Rot_foto2(1,2) + Ny_foto2*Rot_foto2(1,3));
    A(4,2) = (c/D_foto2^2)*(-D_foto2*Rot_foto2(2,2) + Ny_foto2*Rot_foto2(2,3));
    A(4,3) = (c/D_foto2^2)*(-D_foto2*Rot_foto2(3,2) + Ny_foto2*Rot_foto2(3,3));
   
    %Vector de fecho
    w(1,1) = (x0 - c*(Nx_foto1/D_foto1)) - l(1);   
    w(2,1) = (y0 - c*(Ny_foto1/D_foto1)) - l(2);
    w(3,1) = (x0 - c*(Nx_foto2/D_foto2)) - l(3);
    w(4,1) = (y0 - c*(Ny_foto2/D_foto2)) - l(4);
    
    %Vector de correcções aos parâmetros desconhecidos
    delta = -inv(A'*Pl*A)*A'*Pl*w;
    
    %Vector de parâmetros desconhecidos ajustados
    resultado = resultado + delta;
    
    %Vector de Residuos
    v = A*delta + w;
    
    %Vector de observações ajustadas
    l = l + v;
    
    interaccao = interaccao + 1;
end
%----------------------------FIM DO 1º Ciclo-------------------------------

%Matriz de variâncias/covariâncias dos parâmetros desconhecidos
C_delta = sigma0*inv(A'*Pl*A);
C_x = C_delta;

%Matriz de variâncias/covariâncias dos resíduos
C_v = sigma0*(inv(Pl) - A*inv(A'*Pl*A)*A');

%Matriz de variâncias/covariâncias das observações
C_l = Cl-C_v;

%% Escrita dos valores ajustados em ficheiro
if (id_sinal1(k) == 100)
    fprintf(resultados,'Curva perigosa a direita\t\t\t');
    for l = 1:size(resultado,1) 
        fprintf(resultados,'%2.3f\t\t',resultado(l));
    end
    fprintf(resultados,'\n');
    elseif(id_sinal1(k) == 101)
        fprintf(resultados,'Via com entroncamento\t\t\t\t');
        for l = 1:size(resultado,1) 
            fprintf(resultados,'%2.3f\t\t',resultado(l));
        end
        fprintf(resultados,'\n');
    elseif(id_sinal1(k) == 102)
        fprintf(resultados,'Cedencia de prioridade\t\t\t');
        for l = 1:size(resultado,1) 
            fprintf(resultados,'%2.3f\t\t',resultado(l));
        end
        fprintf(resultados,'\n');
    elseif(id_sinal1(k) == 103)
        fprintf(resultados,'Obrigatorio seguir em frente\t\t');
        for l = 1:size(resultado,1) 
            fprintf(resultados,'%2.3f\t\t',resultado(l));
        end
        fprintf(resultados,'\n');
    elseif(id_sinal1(k) == 104)
        fprintf(resultados,'Proibido virar a esquerda\t\t\t');
        for l = 1:size(resultado,1) 
            fprintf(resultados,'%2.3f\t\t',resultado(l));
        end
        fprintf(resultados,'\n');
    elseif(id_sinal1(k) == 105)
        fprintf(resultados,'Proibido ultrapassar\t\t\t\t');
        for l = 1:size(resultado,1) 
            fprintf(resultados,'%2.3f\t\t',resultado(l));
        end
        fprintf(resultados,'\n');
    elseif(id_sinal1(k) == 106)
        fprintf(resultados,'Proibido circular a >40km/h\t\t\t');
        for l = 1:size(resultado,1) 
            fprintf(resultados,'%2.3f\t\t',resultado(l));
        end
        fprintf(resultados,'\n');
    elseif(id_sinal1(k) == 107)
        fprintf(resultados,'Proibido circular a >60km/h\t\t\t');
        for l = 1:size(resultado,1) 
            fprintf(resultados,'%2.3f\t\t',resultado(l));
        end
        fprintf(resultados,'\n');
    elseif(id_sinal1(k) == 108)
        fprintf(resultados,'Proibido circular a >80km/h\t\t\t');
        for l = 1:size(resultado,1) 
            fprintf(resultados,'%2.3f\t\t',resultado(l));
        end
        fprintf(resultados,'\n');
    end
end
%% Fecho de ficheiros
fclose(resultados);