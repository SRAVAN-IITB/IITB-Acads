library ieee;
use ieee.std_logic_1164.all;
entity ALU_MAX_AND_NOT_Eq is
port (
A: in std_logic_vector(3 downto 0);
B: in std_logic_vector(3 downto 0);
sel: in std_logic_vector(1 downto 0);
op: out std_logic_vector(7 downto 0)
) ;
end ALU_MAX_AND_NOT_Eq;

architecture a1 of ALU_MAX_AND_NOT_Eq is

function sub(A: std_logic_vector	:= "0000"; B: std_logic_vector := "0000") 
				return std_logic_vector is
		variable carry: std_logic_vector(8 downto 0);
		variable B_new: std_logic_vector(7 downto 0);
		variable Diff: std_logic_vector(7 downto 0);
		variable A_new: std_logic_vector(7 downto 0);
		variable B_temp: std_logic_vector(7 downto 0);
		begin
			A_new := "0000" & A;
			B_temp := "0000" & B;
			for i in 0 to 7 loop
				B_new(i) := not (B_temp(i));
			end loop;
				
			carry(0) := '1';

			for i in 0 to 7 loop
					Diff(i) := A_new(i) xor B_new(i) xor carry(i);
					carry(i+1) := (((A_new(i) xor B_new(i)) and carry(i)) or (A_new(i) and B_new(i)));
		 	end loop;
	return Diff;
end function sub;

function isequal(A: in std_logic_vector(3 downto 0) :="0000"; 
						B: in std_logic_vector(3 downto 0):="0000")
	return std_logic_vector is
	variable AxnorB: std_logic_vector(3 downto 0):=(others =>'0');
	begin
		for i in 0 to 3 loop
				AxnorB(i) := A(i) xnor B(i);
		end loop;
		if AxnorB = "1111"
		then 
			return "00" & A;
		else
			return "000000";
		end if;
	end function isequal;

function MAX(A: in std_logic_vector(3 downto 0) :="0000")
	return std_logic_vector is
	variable Diff: std_logic := '0';
	variable Eq: std_logic_vector(5 downto 0) := "000000"
	begin
		Eq := isequal(A,B);
		Diff := sub(A,B);
			if E /= "111111"
			then 
				return "000000";
			else
				if Diff = '1'
				then
					return ("00" & A);
				else
					return ("00" & B);
				end if;
			end if;
	end function MAX;




begin
alu : process( A, B, sel )
begin
		case sel is
			when "00" =>
				op <= MAX(A,B);
			when "01" =>
				op <= ANDing(A,B);
			when "10" =>
				op <= compl(A);
			when others =>
				op <= isequal(A,B);
		end case;
	end process ; --alu
end a1 ; -- a1