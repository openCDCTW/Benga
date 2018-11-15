import React from 'react';
import ReactDOM from 'react-dom';

export default class Options extends React.Component {


    render(){
        return (
        	<form name="scheme">
                <br />
                <div>
                <label>Database:&nbsp;&nbsp;</label>
                    <select name="database">
                    <option value="Salmonella_enterica">Salmonella_enterica</option>
                    <option value="Vibrio_cholerae">Vibrio cholerae</option>
                    <option value="Listeria_monocytogenes">Listeria monocytogenes</option>
                    <option value="Campylobacter_jejuni">Campylobacter jejuni</option>
                    </select>
                </div>
                <br />
                <div>
                <label>Email:&nbsp;&nbsp;&nbsp; </label>
                <input name="email" type="text" placeholder="Input your Email here"/>
                </div>
        	</form>
        );
    }
}