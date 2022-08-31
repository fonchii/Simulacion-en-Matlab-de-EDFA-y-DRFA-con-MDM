function modos = modesparam(modes,wavelength,fibra)
% fibra.AN = 0.185;   fibra.AN = 0.2 ; fibra.n1 = 1.45 ; fibra.n2 = sqrt((fibra.n1^2-fibra.AN^2)); fibra.radio = 15e-6
% signal.modos = ['01'] ; 
% modos = modesparam(signal.modos,1550e-9,fibra)

    coreRadius   = fibra.radio*1e6 ;
    nCore        = fibra.n1 ;
    nCladding    = fibra.n2 ;
    lambdas      = wavelength.*1e6 ;
    modos = struct ;

    [solution, totalNumModes] = find_LP_modes(coreRadius, nCore, nCladding, lambdas) ; 
    %for i = 1:length(modes) % modes es un solo str!
    for j = 1:length(solution)
        %mmi1 = extract( modes(i) ,1) ; mmj1 = extract( modes(i) ,1);
        %if solution(j).l == str2num( extract( modes(i) ,1) ) && solution(j).m == str2num( extract( modes(i) ,2) )
        if solution(j).l == str2num( modes(1) ) && solution(j).m == str2num( modes(2) )
            modos.modo = modes ;% modos(i).modo = modes(i) ;
            modos.u = solution(j).u; %modos(i).u = solution(j).u ;
            modos.w = solution(j).w ; %modos(i).w = solution(j).w ;
            modos.beta = solution(j).beta ;%modos(i).beta = solution(j).beta ;
            modos.allsolution = solution ;%modos(i).allsolution = solution ;
        end
    end
    %end
    
    %if length(modes) ~= length(modos) || isfield(modos, 'modo') == 0
    if isfield(modos, 'modo') == 0
        modos = -1;
        error("No es posible transmitir el modo en la fibra" )
    end
    