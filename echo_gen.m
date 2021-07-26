function echo_sound = echo_gen(input, fs, delay, amp)    
    new_sr = round(fs*delay);
    no_echo = [input; zeros(new_sr,1)];
    echo_effect = [zeros(new_sr,1); input*amp];
    echo_sound = no_echo + echo_effect;
    
    norm_factor = max(abs(echo_sound));
    if norm_factor > 1
        echo_sound = echo_sound./norm_factor;
    end

   