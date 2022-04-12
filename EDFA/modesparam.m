function modos = modesparam(modes,wavelength,fibra)
%  modos = modesparam(signal.modos,1520e-9,fibra)
    coreRadius   = fibra.radio*1e6 ;
    nCore        = fibra.n1 ;
    nCladding    = fibra.n2 ;
    lambdas      = wavelength.*1e6 ;
    modos = struct ;

    [solution, totalNumModes] = find_LP_modes(coreRadius, nCore, nCladding, lambdas) ; 
    for i = 1:length(modes)
        for j = 1:length(solution)
            if solution(j).l == str2num( extract( modes(i) ,1) ) && solution(j).m == str2num( extract( modes(i) ,2) )
                modos(i).modo = modes(i) ;
                modos(i).u = solution(j).u ;
                modos(i).w = solution(j).w ;
                modos(i).beta = solution(j).beta ;
                modos(i).allsolution = solution ;
            end
        end
    end
    
    if length(modes) ~= length(modos) || isfield(modos, 'modo') == 0
        modos = -1;
        error("No es posible transmitir el modo en la fibra" )
    end
    