import React from 'react';
import ReactDOM from 'react-dom';

export default class Options_Database extends React.Component {

    render(){
        return (
        		<div>
        		<label>Database:&nbsp;&nbsp;</label>
        			<select name="choose_database" id="database">
        			<option value="Salmonella_7k">S.enterica (7k samples)</option>
        			<option value="Vibrio_cholerae">Vibrio cholerae</option>
        			<option value="Listeria_monocytogenes">Listeria monocytogenes</option>
        			<option value="Campylobacter_jejuni">Campylobacter jejuni</option>
        			</select>
        		</div>
        );
    }
}