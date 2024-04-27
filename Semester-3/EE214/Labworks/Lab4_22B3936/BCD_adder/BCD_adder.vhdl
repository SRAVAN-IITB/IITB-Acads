library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity BCD_Adder is
port (A_3, A_2, A_1, A_0, B_3, B_2, B_1, B_0: in std_logic; Coutt, S_3, S_2, S_1, S_0 : out std_logic);
end entity BCD_Adder;

architecture Struct of BCD_Adder is
--TO BE CHECKED!
signal zero,cout,  S_3dash, S_2dash, S_1dash, S_0dash, C_dash, andS32, andS31, donokaOR, Z : std_logic;

component fourbit_adder is
	port (A_3, A_2, A_1, A_0, B_3, B_2, B_1, B_0: in std_logic; Cout, S_3, S_2, S_1, S_0 : out std_logic);
end component;

begin
	zero <= '0';
  --component instances
  
  fourbit_adder_0 : fourbit_adder port map (A_3 => A_3, A_2 => A_2, A_1 => A_1, A_0 => A_0, B_3 => B_3, B_2 => B_2, B_1 => B_1, B_0 => B_0, Cout => C_dash, S_3 => S_3dash, S_2 => S_2dash, S_1 => S_1dash, S_0 => S_0dash);
  fourbit_adder_1 : fourbit_adder port map (A_3 => S_3dash, A_2 => S_2dash, A_1 => S_1dash, A_0 => S_0dash, B_3 => zero, B_2 => Cout, B_1 => Cout, B_0 => zero, S_3 => S_3, S_2 => S_2, S_1 => S_1, S_0 => S_0);
  AND1 : AND_2 port map (A => S_3dash, B => S_2dash, Y => andS32);
  AND2 : AND_2 port map (A => S_3dash, B => S_1dash, Y => andS31);
  OR1 : OR_2 port map (A => andS32, B => andS31, Y => donokaOR);
  OR2 : OR_2 port map (A => donokaOR, B => C_dash, Y => Cout);
  coutt <= cout;
  end Struct;