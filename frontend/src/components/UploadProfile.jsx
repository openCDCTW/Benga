import React from 'react';
import ReactDOM from 'react-dom';
import Navigation from './Navigation.jsx';
import DropzoneComponent from 'react-dropzone-component';
import { Link } from 'react-router-dom';
import Button from '@material-ui/core/Button';
import { withStyle } from '@material-ui/core/styles';
import CloudUploadIcon from '@material-ui/icons/CloudUpload';
import Icon from '@material-ui/core/Icon';
import DeleteIcon from '@material-ui/icons/Delete';

class Upload_profile extends React.Component {

    constructor(props) {

        super(props);

        fetch('api/dendrogram/upload/', {method:'POST'})
        .then(function(res){
           return res.json();
        }).then(function(batch){
           return getID(batch);
        });

        var getID=function(data){
            window.batchid = data.id;
        };

        this.state = {
            to: "/upload_profile",
            upload_confirm: false,
        };

        this.djsConfig = {
            dictDefaultMessage:"Drop files or click to upload files",
            addRemoveLinks: true,
            acceptedFiles: ".tsv",
            autoProcessQueue: false,
            parallelUploads: 200,
            init:function(){
                this.on("sending", function(file, xhr, formData){
                    formData.append("batch_id", window.batchid);
                });

            }
        }

        this.componentConfig = {
            iconFiletypes: ['.tsv'],
            showFiletypeIcon: true,
            postUrl: 'api/dendrogram/profile/'
        };

        this.dropzone = null;
    }

    handlePost() {

        let fileCheck = this.dropzone.files.length;

        if(fileCheck < 5){
            alert('Please upload at least 5 files');
            return ;
        }

        this.dropzone.processQueue();
        this.setState(state => ({ to: '/dendrogram_view', upload_confirm: true}));
    
    }

    submit(){

        if (this.state.upload_confirm == false ){
            alert('Please upload files first! (At least 5 files)');
            return ;
        }

        var scheme = {};
        scheme.id = window.batchid;
        fetch('api/dendrogram/plot/', {
            method:'POST',
            headers: new Headers({'content-type': 'application/json'}),
            body: JSON.stringify(scheme)
        });
    }

    remove(){

        this.dropzone.removeAllFiles();
        
        this.setState(state => ({ to: '/upload_profile', upload_confirm: false}));
        fetch('api/dendrogram/upload/', {method:'POST'})
        .then(function(res){
           return res.json();
        }).then(function(batch){
           return getID(batch);
        });

        var getID=function(data){
            window.batchid = data.id;
        };
    }

    render() {
        const config = this.componentConfig;
        const djsConfig = this.djsConfig;
        const eventHandlers = {
            init: dz => this.dropzone = dz,
        }

        return (
            <div>
                <Navigation value={1}/>
                <br />
                <br />
                <DropzoneComponent config={config} eventHandlers={eventHandlers} djsConfig={djsConfig} />
                <br />
                <br />
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <Button variant="contained" color="default" onClick={this.handlePost.bind(this)}>
                        Upload
                        &nbsp;&nbsp;
                        <CloudUploadIcon />
                    </Button>
                    &nbsp;&nbsp;&nbsp;&nbsp;
                    <Button variant="contained" color="secondary" onClick={this.remove.bind(this)}>
                        Remove all files
                        &nbsp;&nbsp;
                        <DeleteIcon />
                    </Button>
                    &nbsp;&nbsp;&nbsp;&nbsp;
                    <Link to={this.state.to} style={{ textDecoration:'none' }}>
                        <Button variant="contained" color="primary" onClick={this.submit.bind(this)}>
                            Draw &nbsp; dendrogram
                            &nbsp;&nbsp;
                            <Icon>send</Icon>
                        </Button>
                    </Link>
                </div>
                <br />
                <br />
                <br />
            </div>
        );
    }
}

export default Upload_profile;