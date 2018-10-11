import React from 'react';
import ReactDOM from 'react-dom';
import Scheme from './Options_Scheme.jsx';
import Database from './Options_Database.jsx';
import Email from './Options_email.jsx';

export default class Options extends React.Component {

    // constructor(){
    //     super();
    //     this.handleSubmit = this.handleSubmit.bind(this);
    // }

    // handleSubmit(event){
    //     event.preventDefault();
    //     // const data = new FormData(event.target);
    //     let id = '';
    //     // var data={}
        

    //     fetch('api/profiling/upload/',{
    //         method:'POST',
    //     }).then(function(res){
    //        return res.json();
    //     }).then(function(create_uuid){
    //         var data={
    //             'batch_id': create_uuid.id,
    //             'filename':'test',
                
    //         };
    //         return data;
    //     }).then(function(data){

    //         fetch('api/profiling/sequence/',{
    //         method:'POST',
    //         headers:{
    //             "Content-Type":"application/json",
    //             },
    //         body:JSON.stringify(data),
    //         })
    //     });


    // }


    render(){
        return (
        	<form onSubmit={this.handleSubmit}>
        		<div className="Scheme">
        			<Scheme />
                </div>
                <br />
                <div className="Database">
                	<Database />
                </div>
                <br />
                <div className="Email">
                	<Email />
                </div>
                <br />
                <button> send </button>
        	</form>
        );
    }
}