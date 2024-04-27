library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity Fourbit_Adder is
port (A_3, A_2, A_1, A_0, B_3, B_2, B_1, B_0: in std_logic; Cout, S_3, S_2, S_1, S_0 : out std_logic);
end entity Fourbit_Adder;

architecture Struct of Fourbit_Adder is
--TO BE CHECKED!
signal Cout_0, Cout_1, Cout_2 : std_logic;

component Full_adder is
port (A, B, Cin: in std_logic; S, Cout: out std_logic);
end component;

begin
  --component instances
  
  Half_adder_0 : HALF_ADDER port map (A => A_0, B => B_0, S => S_0, C => Cout_0);
  Full_adder_1 : Full_adder port map (A => A_1, B => B_1, Cin => Cout_0, S => S_1, Cout => Cout_1);
  Full_adder_2 : Full_adder port map (A => A_2, B => B_2, Cin => Cout_1, S => S_2, Cout => Cout_2);
  Full_adder_3 : Full_adder port map (A => A_3, B => B_3, Cin => Cout_2, S => S_3, Cout => Cout);
  
  end Struct;