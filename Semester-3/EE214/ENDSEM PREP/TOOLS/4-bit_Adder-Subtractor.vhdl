library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity Fourbit_AdderSubtractor is
port (A3, A2, A1, A0, B3, B2, B1, B0, M: in std_logic; Cout, S3, S2, S1, S0 : out std_logic);
end entity Fourbit_AdderSubtractor;

architecture Struct of Fourbit_AdderSubtractor is
--TO BE CHECKED!
signal XOR0, XOR1, XOR2, XOR3, Cout_0, Cout_1, Cout_2 : std_logic;

component Full_adder is
port (A, B, C: in std_logic; Sum, Cout: out std_logic);
end component;

begin
  --component instances
  XORR1 : XOR_2 port map (A => M, B => B0, Y => XOR0);
  XORR2 : XOR_2 port map (A => M, B => B1, Y => XOR1);
  XORR3 : XOR_2 port map (A => M, B => B2, Y => XOR2);
  XORR4 : XOR_2 port map (A => M, B => B3, Y => XOR3);
  Full_adder_0 : Full_adder port map (A => XOR0, B => A0, C => M, Sum => S0, Cout => Cout_0);
  Full_adder_1 : Full_adder port map (A => XOR1, B => A1, C => Cout_0, Sum => S1, Cout => Cout_1);
  Full_adder_2 : Full_adder port map (A => XOR2, B => A2, C => Cout_1, Sum => S2, Cout => Cout_2);
  Full_adder_3 : Full_adder port map (A => XOR3, B => A3, C => Cout_2, Sum => S3, Cout => Cout);
  
  end Struct;