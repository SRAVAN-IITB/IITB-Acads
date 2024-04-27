library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity Full_adder is
port (A, B, C: in std_logic; Sum, Cout: out std_logic);
end entity Full_adder;

architecture Struct of Full_adder is
signal AxorB, AandB, BandCin, AandCin, orr : std_logic;

begin
  --component instances
  XOR1 : XOR_2 port map (A => A, B => B, Y => AxorB);
  XOR2 : XOR_2 port map (A => AxorB, B => C, Y => Sum);

  AND1 : AND_2 port map (A => A, B => B, Y => AandB);
  AND2 : AND_2 port map (A=>C, B=>B, Y=>BandCin);
  AND3 : AND_2 port map (A=>C, B=>A, Y=>AandCin);
  OR1 : OR_2 port map (A=>AandB, B=>BandCin, Y=>orr);
  OR2 : OR_2 port map (A=>orr, B=>AandCin, Y=>Cout);

  end Struct;
