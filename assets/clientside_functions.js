window.dash_clientside = Object.assign({}, window.dash_clientside, {
    login: {
        password_strength_zxcvbn: function (password) { // output: value, color
            let result = zxcvbn(password)
            let guesses_log10 = result['guesses_log10']
            let value = Math.min(100, ( guesses_log10 / 16 )*100)
            if(value < 40){
                return [value, 'red']
            }else if(value < 60){
                return [value, 'orange']
            }else if(value < 80){
                return [value, 'lime']
            }else{
                return [value, 'green']
            }

        }
    },

    overview: {
        grid_container_resize: function(width, height){
            window.dispatchEvent(new Event('resize'))
            return height
        },

        update_PlotPanel_display_idx: function (list_cur_uuid){
            // console.log(window.dash_clientside.callback_context)
            // let out_uuid = window.dash_clientside.callback_context.outputs_list.map(
            //     x => x['id']['index']
            // )
            let out_idx = Array.from(list_cur_uuid.keys())
            let out_label = out_idx.map(x => 'Panel '+(x+1))
            return out_label
        }
    }
});
