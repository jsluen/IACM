function f=normalize01(f)

    fmin  = min(f(:));
    fmax  = max(f(:));
    f = (f-fmin)/(fmax-fmin); 
    
end
