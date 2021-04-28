function ss=add_swarm_to_XDMF(field_name,file_name,is_vector,dimensions,is_int)


nl = sprintf('\n'); % new line
if is_vector
    d='3';
    type_vec_scal='Vector';
else
    d='1';
    type_vec_scal='Scalar';
end

if is_int
    type='"Int"';
else
    type='"Float" Precision="8"';
end

ss=['   <Attribute Type="',type_vec_scal,'" Center="Node" Name="',field_name,'">',nl, ...
'           <DataItem Format="HDF" NumberType=',type,' Dimensions="',dimensions,32,d,'">',file_name,'</DataItem>',nl, ...
'   </Attribute>',nl];

% '         <Attribute Type="Scalar" Center="Node" Name="materialSwarm-MaterialIndex">',nl,...
% '            <DataItem Format="HDF" NumberType="Int" Dimensions="180245 1">../materialSwarm.00187.1of64.h5:/materialSwarm-MaterialIndex</DataItem>',nl,...
% '         </Attribute>',nl,nl,...
end