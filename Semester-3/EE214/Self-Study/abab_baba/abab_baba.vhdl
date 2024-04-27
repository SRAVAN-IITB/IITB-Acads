
library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Flipflops.all;

entity seqDetector is
port(reset, inp, clock : in std_logic; --inp=0 => a inp=1 => b
        y : out std_logic);
end entity seqDetector;

architecture struct of seqDetector is

    signal d0, d1, d2, q0, q1, q2, q0_bar, q1_bar, q2_bar, inp_bar : std_logic;
begin
    q0_bar <= not q0;
    q1_bar <= not q1;
    q2_bar <= not q2;
    inp_bar <= not inp;
    
    
    ------------------Here we are defining the combinational block L1 which is the input to the FF0----------------
    d0 <= (not reset) and ((inp_bar and q1_bar) or (q2_bar and q1_bar and q0) 
    or (q2_bar and q0 and inp_bar) or (q2 and q0_bar and inp_bar));
    

    ------------------D Flip Flop 0--------------------------------
    dFlipFlop_0 : dff_reset port map(d0, clock, reset, q0);
    
    
    
    ------------------Here we are defining the combinational block L2 whoich is the input to the FF1----------------
    d1 <= (not reset) and ((q2_bar and inp) or (q1_bar and inp));
    
    ------------------D Flip Flop 1--------------------------------
    dFlipFlop_1 : dff_reset port map(d1, clock, reset, q1);
    
    
    
    ------------------Here we are defining the combinational block L3 whoich is the input to the FF2----------------
    d2 <= (not reset) and ((q2_bar and inp_bar and q1) or (q0_bar and inp_bar and q1) or (q2 and q1_bar and inp));
    
    ------------------D Flip Flop 2--------------------------------
    dFlipFlop_2 : dff_reset port map(d2, clock, reset, q2);


    ------------------Here goes the output------------------------------
    y<= q2 and ((q1_bar and q0 and inp) or (q1 and q0_bar and inp_bar));

end struct;
