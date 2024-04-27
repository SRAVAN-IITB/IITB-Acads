library ieee;
use ieee.std_logic_1164.all;
library work;
use work.Gates.all;

entity mult is
port (a: in std_logic_vector( 3 downto 0);
		b: in std_logic_vector( 2 downto 0); 
		o: out std_logic_vector(6 downto 0));
end entity mult;

architecture Struct of mult is
	component Fourbit_Adder is
	port (A_3, A_2, A_1, A_0, B_3, B_2, B_1, B_0: in std_logic; 
			Cout, S_3, S_2, S_1, S_0 : out std_logic);
	end component;
	
	signal s: std_logic_vector (11 downto 0);
	signal x: std_logic_vector (3 downto 0);
begin
  --component instances
  
  i_loop : for i in 0 to 3 generate
  j_loop : for j in 0 to 2 generate
  ANDing: AND_2 port map (A=> a(i), B=>b(j), Y=>s(i + 4*j) );
  end generate j_loop;
  end generate i_loop;
  
  -- x(3) <= '0';
  fouradd1: fourbit_Adder port map (A_3 => '0', A_2 => s(3), A_1 => s(2), A_0 => s(1), 
												B_3 => s(7), B_2 => s(6), B_1 => s(5), B_0 => s(4), 
												Cout => x(3), S_3 => x(2), S_2 => x(1), S_1 => x(0), S_0 => o(1));

  fouradd2: fourbit_Adder port map (A_3 => x(3), A_2 => x(2), A_1 => x(1), A_0 => x(0), 
												B_3 => s(11), B_2 => s(10), B_1 => s(9), B_0 => s(8), 
												Cout => o(6), S_3 => o(5), S_2 => o(4), S_1 => o(3), S_0 => o(2));
  o(0) <= s(0); 
  end Struct;