library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity Half_subtractor is
port (A, B: in std_logic; Dif, Bor: out std_logic);
end entity Half_subtractor;

architecture Struct of Half_subtractor is
signal mid1, mid2 : std_logic;

component XOR_gate is
port (A, B: in std_logic; Y: out std_logic);
end component;

begin
  --component instances
  XOR1 : XOR_gate port map (A => A, B => B, Y => Dif);
  NOR1 : NOR_2 port map (A => B, B => B, Y => mid1);
  NOR2 : NOR_2 port map (A=>mid1, B=>A, Y=>Bor);
  end Struct;