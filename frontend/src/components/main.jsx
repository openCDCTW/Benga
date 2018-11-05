import React from 'react';
import ReactDOM from 'react-dom';
import DropzoneComponent from 'react-dropzone-component';
import Header from './header.jsx';
import Footer from './footer.jsx';
import Options from './Options.jsx';

// Render

class Main extends React.Component {

    constructor(props) {

        super(props);

        fetch('api/profiling/upload/',{method:'POST'})
        .then(function(res){
           return res.json();
        }).then(function(batch){
           return getID(batch);
        });

        var getID=function(data){
            window.batchid=data.id;
        };


        // TODO: poor performnce issue with everytime load the component will
        // fetch API once.

        this.djsConfig = {
            dictDefaultMessage:"Drop files or click to upload files",
            addRemoveLinks: true,
            acceptedFiles: ".fasta,.fa,.fna",
            autoProcessQueue: false,
            parallelUploads:200,
            init:function(){
                this.on("sending",function(file,xhr,formData){
                    formData.append("batch_id",window.batchid);
                });

                this.on("queuecomplete",function(){

                    var scheme = {};
                    scheme.occurrence = document.scheme.occurrence.value;
                    scheme.database = document.scheme.database.value;
                    scheme.id = window.batchid;

                    fetch('api/profiling/profiling/',{
                        method:'POST',
                        headers:new Headers({'content-type':'application/json'}),
                        body:JSON.stringify(scheme)
                    });
                });
            }
        }

        this.componentConfig = {
            iconFiletypes: ['.fasta','.fna','.fa'],
            showFiletypeIcon: true,
            postUrl: 'api/profiling/sequence/'
        };

        this.dropzone = null;
    }

    handlePost() {
        this.dropzone.processQueue();
    }

    reload(){
        window.location.reload();
    }

    
    render() {

        const config = this.componentConfig;
        const djsConfig = this.djsConfig;

        // For a list of all possible events (there are many), see README.md!
        const eventHandlers = {
            init: dz => this.dropzone = dz,
        }

        return (
            <div>
                <div>
                <Header />
                </div>
                <br />
                <br />
                <DropzoneComponent config={config} eventHandlers={eventHandlers} djsConfig={djsConfig} />
                <br />
                <Options />
                <br />
                <button onClick={this.handlePost.bind(this)}> Upload </button> &nbsp;&nbsp;&nbsp;&nbsp;
                <button onClick={this.reload}> Reset </button>
                <br />
                <br />
                <div>
                <Footer />
                </div>


            </div>
        );
    }
}

ReactDOM.render(<Main /> ,document.getElementById('host'));
