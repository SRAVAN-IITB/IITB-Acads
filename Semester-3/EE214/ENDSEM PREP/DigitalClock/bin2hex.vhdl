library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity bin2hex is
	port (
		 Bin: in std_logic_vector(3 downto 0);
		 Hout: out std_logic_vector(6 downto 0)
	);
end bin2hex;

architecture Behavioral of bin2hex is
begin
 process(Bin)
 begin
  case(Bin) is
   when "0000" =>  Hout <= "1000000"; --0--
   when "0001" =>  Hout <= "1111001"; --1--
   when "0010" =>  Hout <= "0100100"; --2--
   when "0011" =>  Hout <= "0110000"; --3--
   when "0100" =>  Hout <= "0011001"; --4-- 
   when "0101" =>  Hout <= "0010010"; --5--    
   when "0110" =>  Hout <= "0000010"; --6--
   when "0111" =>  Hout <= "1111000"; --7--   
   when "1000" =>  Hout <= "0000000"; --8--
   when "1001" =>  Hout <= "0010000"; --9--
   when "1010" =>  Hout <= "0001000"; --a--
   when "1011" =>  Hout <= "0000011"; --b--
   when "1100" =>  Hout <= "1000110"; --c--
   when "1101" =>  Hout <= "0100001"; --d--
   when "1110" =>  Hout <= "0000110"; --e--
   when others =>  Hout <= "0001110"; 
   end case;
 end process;
end Behavioral;