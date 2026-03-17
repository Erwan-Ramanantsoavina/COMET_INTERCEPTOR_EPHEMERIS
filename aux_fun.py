def ctr(name , date , x0_loc) :
    ctr = {
        'type': 'control',
        'epoch': str(date),
        'name': name,
        'state': [
            {'name': 'SC_center',
            'point': 'Earth',
            'coi': 'Earth',
            'axes': 'SEROT',
            'project': False,
            'dynamics': 'EMS_gravity',
            'value': {
                'pos_x': str(x0_loc[0]) + ' km',
                'pos_y': str(x0_loc[1]) + ' km',
                'pos_z': str(x0_loc[2]) + ' km',
                'vel_x': str(x0_loc[3]) + ' km/s',
                'vel_y': str(x0_loc[4]) + ' km/s',
                'vel_z': str(x0_loc[5]) + ' km/s'
            }},
            {'name': 'SC_dv',
             'value': '0 m/s'}
        ]}
    return ctr

def man(name , reference , dt_s , dv_vect_kms) :
    man = {'type': 'manoeuvre',
            'name': name ,
            'model': 'impulsive',
            'input': 'SC',
            'thruster' : 'main',
            'config': { 'point': {'reference': reference, 'dt' : str(dt_s) + ' s'} ,
                'direction': {
                           'axes': 'SEROT',
                           'dv_x': str(dv_vect_kms[0]) +' km/s' ,
                           'dv_y': str(dv_vect_kms[1]) +' km/s' ,
                           'dv_z': str(dv_vect_kms[2]) +' km/s'
                       }}}
    return man

def match(name , ref_left , dt_left , ref_right , dt_right) :
    match = {'type': 'match' ,
             'name': name ,
             'input': 'SC' ,
             'left': {
                 'reference': ref_left,
                 'dt': str(dt_left) + ' s'},
             'right': {
                 'reference': ref_right,
                 'dt': str(dt_right) + ' s'},
             'body': 'Earth',
             'vars': 'cart'}
    return match
    