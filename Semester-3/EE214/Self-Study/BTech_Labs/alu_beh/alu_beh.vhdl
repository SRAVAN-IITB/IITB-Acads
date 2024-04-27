library ieee;
use ieee.std_logic_1164.all;
entity alu_beh is
port (
A: in std_logic_vector(3 downto 0);
B: in std_logic_vector(3 downto 0);
sel: in std_logic_vector(1 downto 0);
op: out std_logic_vector(7 downto 0)
) ;
end alu_beh;

architecture a1 of alu_beh is

function full_adder(A, B, C: in std_logic) return std_logic_vector is
variable s: std_logic_vector(1 downto 0);
begin
s(1):= (a or b) and ((b or c) and (c or a));
s(0):= c xor (a xor b);
return s;
end full_adder;

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

function xoring(A: in std_logic_vector(3 downto 0) :="0000"; 
						B: in std_logic_vector(3 downto 0):="0000")
	return std_logic_vector is
	variable AxorB: std_logic_vector(7 downto 0):=(others =>'0');
	begin
		for i in 0 to 3 loop
				AxorB(i) := A(i) xor B(i);
		end loop;
		AxorB(7 downto 4) := "0000";
		return AxorB;
	end function xoring;

function doubleA(A: in std_logic_vector(3 downto 0) :="0000")
	return std_logic_vector is
	variable ATwice: std_logic_vector(7 downto 0):=(others =>'0');
	begin
		for i in 1 to 4 loop
				ATwice(i) := A(i-1);
		end loop;
		ATwice(7 downto 5) := "000";
		ATwice(0) := '0';
		return ATwice;
	end function doubleA;


begin
alu : process( A, B, sel )
begin
	case sel is
		when "00" =>
			op <= A & B;
		when "01" =>
			op <= sub(A,B);
		when "10" =>
			op <= xoring(A,B);
		when others =>
			op <= doubleA(A);
	end case;		
end process ; --alu
end a1 ; -- a1