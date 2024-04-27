library ieee;
use ieee.std_logic_1164.all;
entity ALU_MAX_AND_NOT_Eq is
port (
A: in std_logic_vector(3 downto 0);
B: in std_logic_vector(3 downto 0);
--sel: in std_logic_vector(1 downto 0);
op: out std_logic_vector(5 downto 0)
) ;
end ALU_MAX_AND_NOT_Eq;

architecture a1 of ALU_MAX_AND_NOT_Eq is
signal sel: std_logic_vector(1 downto 0);
	
function isequal(A: in std_logic_vector(3 downto 0) :="0000"; 
						B: in std_logic_vector(3 downto 0):="0000")
	return std_logic_vector is
	variable AxnorB: std_logic_vector(3 downto 0):=(others =>'0');
	begin

		if A = B
		then
			return "00" & A;
		else
			return "000000";
		end if;
	end function isequal;

function MAX(A: in std_logic_vector(3 downto 0) :="0000";
				B: in std_logic_vector(3 downto 0):="0000")
	return std_logic_vector is
	variable Eq: std_logic_vector(5 downto 0):=(others =>'0');

	
	begin
		Eq := isequal(A,B);
			if Eq /= "000000"
			then 
				return "000000";
			else
				if A > B
				then
					return ("00" & A);
				else
					return ("00" & B);
				end if;
			end if;
	end function MAX;

function compl(A: in std_logic_vector(3 downto 0) :="0000"; 
						B: in std_logic_vector(3 downto 0):="0000")
	return std_logic_vector is
	variable Abar: std_logic_vector(3 downto 0):=(others =>'0');
	begin
		for i in 0 to 3 loop
				Abar(i) := not A(i);
		end loop;
			return "00" & Abar;
	end function compl;
	
function ANDing(A: in std_logic_vector(3 downto 0) :="0000"; 
						B: in std_logic_vector(3 downto 0):="0000")
	return std_logic_vector is
	variable AandB: std_logic_vector(3 downto 0):=(others =>'0');
	begin

		for i in 0 to 3 loop
				AandB(i) := A(i) and B(i);
		end loop;
			return "00" & AandB;
	end function ANDing;

begin
	sel <= B(3) & A(3);
alu : process( A, B,sel )
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