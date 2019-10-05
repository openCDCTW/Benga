import React from 'react';
import ReactDOM from 'react-dom';
import DropzoneComponent from 'react-dropzone-component';
import { Link } from 'react-router-dom';
import Options from './Options.jsx';
import Button from '@material-ui/core/Button';
import { withStyles } from '@material-ui/core/styles';
import CloudUploadIcon from '@material-ui/icons/CloudUpload';
import Icon from '@material-ui/core/Icon';
import DeleteIcon from '@material-ui/icons/Delete';
import ReplyIcon from '@material-ui/icons/Reply';
import DownloadIcon from '@material-ui/icons/CloudDownload';
import SearchBar from './SearchBar.jsx';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import blue from '@material-ui/core/colors/blue';
// import { Scrollbars } from 'react-custom-scrollbars';
import download from 'downloadjs';


const styles = theme => ({
    cssRoot:{
        color: theme.palette.getContrastText(blue[900]),
        backgroundColor: blue[900],
        '&:hover': {
            backgroundColor:blue[600],
        },
    }
})

class Upload_contigs extends React.Component {

    constructor(props) {

        super(props);

        window.databaseName = "Vibrio_cholerae";
        window.fileName = [];

        this.state = {
            switch: false,
        };

        this.djsConfig = {
            dictDefaultMessage:"Drag contig file(s) here (up to 100 files)",
            dictRemoveFile:"Remove",
            addRemoveLinks: true,
            maxFilesize: 10,
            acceptedFiles: ".fasta,.fa,.fna",
            autoProcessQueue: false,
            timeout: 0,
            parallelUploads: 100,
            init:function(){
                this.on("addedfile", function(file){
                    if(file.size < 100000){
                        this.removeFile(file);
                    }
                });
                this.on("sending", function(file, xhr, formData){
                    formData.append("batch_id", window.batchid);
                    formData.append("database", window.databaseName);
                    formData.append("occurrence", 95);
                    window.fileName.push(file.name);
                });
                this.on("success", function(file){
                    file._removeLink.remove();
                    delete file._removeLink;
                });
            }
        }

        this.componentConfig = {
            iconFiletypes: ['.fasta','.fna','.fa'],
            showFiletypeIcon: true,
            postUrl: '/cgMLST/api/profiling/sequence/'
        };

        this.dropzone = null;

    }

    handlePost() {
        
        var fileCheck = this.dropzone.files.length;

        if(fileCheck < 1){
            alert('Please upload contigs file first');
            return ;
        }else if(fileCheck > 100){
            alert('Cannot upload more than 100 files');
            return ;
        }

        if(window.databaseName == ""){
            alert('Please choose a database !');
            return ;
        }

        this.setState(state => ({ switch: true }));

        fetch('/cgMLST/api/profiling/upload/', {method:'POST'})
            .then(function(res){
               return res.json();
            }).then(batch => window.batchid = batch.id);

        function submit(){

            if( window.batchid != undefined ){
                this.dropzone.processQueue();

                var scheme = {};
                scheme.seq_num = fileCheck;
                fetch('/cgMLST/api/profiling/upload/' + window.batchid + '/', { 
                    method:'PATCH',
                    headers: new Headers({'content-type': 'application/json'}),
                    body: JSON.stringify(scheme)
                });

                this.props.history.push("/cgMLST/non-release/profile_result");
                clearInterval(interval);
            }
        };

        let interval = setInterval(submit.bind(this),500);
    }

    remove(){
        this.dropzone.removeAllFiles();
        this.setState(state => ({ switch: false }));
        window.fileName.length = 0;
    }

    query(){

        fetch('/cgMLST/api/profiling/profile/' + this.state.querybyID + '/', { method:'GET'})
        .then(function(response){
            if(response.status != 404){
                return response.json();
            }else{
                return response.status;
            }
        }).then(res => this.setState(state => ({ tmp: res })));

        function result(){
            if(this.state.tmp.zip != undefined){
                this.setState(state => ({ profile_result_zip : '/cgMLST/'+this.state.tmp.zip }))
            }else{
                alert("Data not found. Please input correct ID or try again later.");
            }
            clearInterval(interval);
        }
        let interval = setInterval(result.bind(this),50);
    }

    back(){
        this.setState(state => ({ profile_result_zip: undefined }));
    }
    
    render() {

        const config = this.componentConfig;
        const djsConfig = this.djsConfig;
        const eventHandlers = {
            init: dz => this.dropzone = dz,}
        const { classes } = this.props;

        if(this.state.profile_result_zip == undefined){
            return (
            <div>
                <br />
                <br />
                <div>
                    <div style={{ float:'left', marginLeft:'10px', marginTop:'7px' }}>
                        <Options switch={this.state.switch} />
                    </div>
                    <div style={{ float:'right', marginTop:'35px', marginRight:'25px' }}>
                        <Button variant="contained" color="secondary" onClick={this.remove.bind(this)}>
                                Remove all files
                                &nbsp;&nbsp;
                                <DeleteIcon />
                        </Button>
                    </div>
                </div>
                <br />
                <br />
                <br />
                <br />
                <br />
                <div style = {{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                    <DropzoneComponent config={config} eventHandlers={eventHandlers} 
                        djsConfig={djsConfig} />
                </div>
                <br />
                <br />
                <div style = {{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                    <Button variant="contained" className ={classes.cssRoot} 
                     onClick={this.handlePost.bind(this)}>
                        Submit
                        &nbsp;&nbsp;
                        <CloudUploadIcon />
                    </Button>
                </div>
                <br />
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <SearchBar
                        onChange = {(value) => this.setState({ querybyID: value })}
                        onRequestSearch={this.query.bind(this)}
                        placeholder = {"Input ID to get result"}
                        style = {{
                            width: '90%',
                            margin: '0 auto',
                        }}
                    />
                </div>
                <br />
                <br />
                <br />
                <br />
                <br />
            </div>
            );
        }else{
            return (
                    <div>
                        <br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <a download href={this.state.profile_result_zip} 
                             style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Download profiles (.zip)
                                &nbsp;&nbsp;
                                <DownloadIcon />
                                </Button>
                            </a>
                        </div>
                        <br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <Button variant="contained" color="default" onClick={this.back.bind(this)}>
                                <ReplyIcon />
                                &nbsp;&nbsp;
                                Back
                            </Button>
                        </div>
                        <br />
                    </div>
                );
        }
        
    }
}


export default withStyles(styles)(Upload_contigs);
