function field = choose_field(all_fields,numbering)

switch numbering
    case 1
        field = all_fields.phi;
    case 2
        field = all_fields.s;
    case 3
        field = all_fields.c;
    case 4
        field = all_fields.e;
end

end