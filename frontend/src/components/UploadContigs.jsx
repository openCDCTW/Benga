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

        fetch('api/profiling/upload/', {method:'POST'})
        .then(function(res){
           return res.json();
        }).then(function(batch){
           return getID(batch);
        });

        var getID=function(data){
            window.batchid = data.id;
        };

        // TODO: poor performnce issue with everytime load the component will
        // fetch API once.

        window.databaseName = "Vibrio_cholerae";
        window.fileName = [];

        this.state = {
            switch: false,
        };

        this.djsConfig = {
            dictDefaultMessage:"Drag contig file(s) here (up to 100 files)",
            dictRemoveFile:"Remove",
            addRemoveLinks: true,
            maxFilesize:10,
            acceptedFiles: ".fasta,.fa,.fna",
            autoProcessQueue: false,
            parallelUploads: 100,
            init:function(){
                this.on("addedfile", function(on_load_header_data){
                    
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
            postUrl: 'api/profiling/sequence/'
        };

        this.dropzone = null;

    }

    handlePost() {
        
        let fileCheck = this.dropzone.files.length;

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

        this.dropzone.processQueue();

        var scheme = {};
        scheme.seq_num = this.dropzone.files.length;
        //database ,file 
        fetch('api/profiling/upload/' + window.batchid + '/', { 
            method:'PATCH',
            headers: new Headers({'content-type': 'application/json'}),
            body: JSON.stringify(scheme)
        });

        this.props.history.push("/profile_view");
    }

    remove(){
        
        this.dropzone.removeAllFiles();
        this.setState(state => ({ switch: false }));
        window.fileName.length = 0;

        fetch('api/profiling/upload/', {method:'POST'})
        .then(function(res){
           return res.json();
        }).then(function(batch){
           return getID(batch);
        });

        var getID=function(data){
            window.batchid = data.id;
        };
    }

    query(){

        fetch('api/profiling/profile/' + this.state.querybyID + '/', { method:'GET'})
        .then(response => response.json())
        .then(result => this.setState(state => ({
                profile_result_zip: result.zip })))
        .catch(error => alert("Not found"));

        // fetch('api/profiling/profile/' + this.state.querybyID + '/', { method:'GET'})
        //     .then(response => response.json())
        //     .then(result => this.setState(state => ({
        //         profile_result_zip: result.zip })));
    }

    upload_example_data(){

        let exampleFile = [
            { name:"V.cholerae_01.fa" },
            { name:"V.cholerae_02.fa" },
            { name:"V.cholerae_03.fa" },
            { name:"V.cholerae_04.fa" },
            { name:"V.cholerae_05.fa" },
            { name:"V.cholerae_06.fa" },
        ];

        var scheme = {};
        scheme.seq_num = 6;
        fetch('api/profiling/upload/' + window.batchid + '/', { 
            method:'PATCH',
            headers: new Headers({'content-type': 'application/json'}),
            body: JSON.stringify(scheme)
        });

        let i = 0;
        for(i; i < exampleFile.length; i++){
            this.dropzone.emit("addedfile", exampleFile[i]);
            this.dropzone.emit("success", exampleFile[i]);
            this.dropzone.emit("complete", exampleFile[i]);
            this.dropzone.files.push(exampleFile[i]);
            window.fileName.push(exampleFile[i].name);
        };

        let encodeExampleData = [
            require('./static/Example_data/V.cholerae_01.fa'), 
            require('./static/Example_data/V.cholerae_02.fa'), 
            require('./static/Example_data/V.cholerae_03.fa'), 
            require('./static/Example_data/V.cholerae_04.fa'), 
            require('./static/Example_data/V.cholerae_05.fa'), 
            require('./static/Example_data/V.cholerae_06.fa'), 
        ];

        let decodeExampleData = [];
        let j = 0;

        for(j; j < encodeExampleData.length; j++){
            let tmp = encodeExampleData[j].substring(13,encodeExampleData[j].length);
            tmp = window.atob(tmp);
            decodeExampleData.push(tmp);
        };

        let VC01 = new File([decodeExampleData[0]],'V.cholerae_01.fa');
        let VC02 = new File([decodeExampleData[1]],'V.cholerae_02.fa');
        let VC03 = new File([decodeExampleData[2]],'V.cholerae_03.fa');
        let VC04 = new File([decodeExampleData[3]],'V.cholerae_04.fa');
        let VC05 = new File([decodeExampleData[4]],'V.cholerae_05.fa');
        let VC06 = new File([decodeExampleData[5]],'V.cholerae_06.fa');

        let decodeExampleFile = [ VC01, VC02, VC03, VC04, VC05, VC06 ];

        let k = 0;
        let upload_status = [];

        window.databaseName = "Vibrio_cholerae";

        for(k; k < decodeExampleFile.length; k++){
            let form = new FormData();
            form.append('file',decodeExampleFile[k]);
            form.append('batch_id',window.batchid);
            form.append('occurrence',"95");
            form.append('database',window.databaseName);

            fetch('api/profiling/sequence/', {
                method:'POST',
                body:form ,
            }).then(res => upload_status.push(res.status));
        };
      
        this.setState(state => ({ switch: true }));

        this.props.history.push("/profile_view");

        // function trigger_celery(){

        //     if( upload_status.length == 12 ){
        //         let scheme = {};
        //         scheme.occurrence = "95";
        //         scheme.database = window.databaseName;
        //         scheme.id = window.batchid;
        //         fetch('api/profiling/profiling-tree/', {
        //             method:'POST',
        //             headers: new Headers({'content-type': 'application/json'}),
        //             body: JSON.stringify(scheme)
        //         });
        //         clearInterval(interval);
        //     }
        // };

        // let interval = setInterval(trigger_celery,1500);
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
                    <Button variant="contained" color="default" 
                     onClick={this.upload_example_data.bind(this)}>
                        Example
                        &nbsp;&nbsp;
                        <CloudUploadIcon />
                    </Button>
                    &nbsp;&nbsp;&nbsp;&nbsp;
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
                        placeholder = {"ID to get result"}
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
                        <br />
                    </div>
                );
        }
        
    }
}


export default withStyles(styles)(Upload_contigs);
