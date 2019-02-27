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
import green from '@material-ui/core/colors/green';
// import { Scrollbars } from 'react-custom-scrollbars';
import download from 'downloadjs';


const styles = theme => ({
    cssRoot:{
        color: theme.palette.getContrastText(green[600]),
        backgroundColor: green[500],
        '&:hover': {
            backgroundColor:green[600],
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

        window.databaseName = "";
        window.fileName = [];

        this.state = {
            to: "/",
            upload_confirm: false,
            switch: false,
        };

        this.djsConfig = {
            dictDefaultMessage:"Drop or click to upload contig files",
            addRemoveLinks: true,
            maxFilesize:10,
            acceptedFiles: ".fasta,.fa,.fna",
            autoProcessQueue: false,
            parallelUploads: 200,
            init:function(){
                this.on("addedfile", function(on_load_header_data){
                    // var fileSizeElement = on_load_header_data.previewElement.querySelector(".dz-size");
                    // fileSizeElement.parentNode.removeChild(fileSizeElement);
                });
                this.on("sending", function(file, xhr, formData){
                    formData.append("batch_id", window.batchid);
                });
                this.on("success", function(file){
                    window.fileName.push(file.name);
                    file._removeLink.remove();
                    delete file._removeLink;
                });
                this.on("removedfile", function(file){
                    // file.previewElement.fadeOut('slow');
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
        }

        if(window.databaseName == ""){
            alert('Please choose a database !');
            return ;
        }

        this.dropzone.processQueue();
        this.setState(state => ({ upload_confirm: true , to: '/profile_view' ,
            switch: true }));
    }

    submit(){

        if (this.state.upload_confirm == false){
            alert('Please upload files first! (At least 5 files)');
            return ;
        }

        var scheme = {};
        scheme.occurrence = "95";
        scheme.database = window.databaseName;
        scheme.id = window.batchid;
        fetch('api/profiling/profiling-tree/', {
            method:'POST',
            headers: new Headers({'content-type': 'application/json'}),
            body: JSON.stringify(scheme)
        });

        window.tabSwitch = true;

    }

    remove(){
        
        this.dropzone.removeAllFiles();
        this.setState(state => ({ to: '/', upload_confirm: false, switch: false }));
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

        fetch('api/profiling/profile/' + this.state.querybyID, { method:'GET'})
        .then(response => response.json())
        .catch(error => alert("Data not found"));

        fetch('api/profiling/profile/' + this.state.querybyID, { method:'GET'})
            .then(response => response.json())
            .then(result => this.setState(state => ({
                profile_result_all: result.file,
                profile_result_zip: result.zip,})));

        fetch('api/dendrogram/dendrogram/' + this.state.querybyID, { method: 'GET'})
            .then(response => response.json())
            .then(result => this.setState(state => ({
                png_file: result.png_file, 
                pdf_file: result.pdf_file,
                svg_file: result.svg_file, 
                emf_file: result.emf_file, 
                newick_file: result.newick_file, })));
    }

    upload_example_data(){

        let exampleFile = [
            { name:"V.cholerae_01.fa" },
            { name:"V.cholerae_02.fa" },
            { name:"V.cholerae_03.fa" },
            { name:"V.cholerae_04.fa" },
            { name:"V.cholerae_05.fa" },
            { name:"V.cholerae_06.fa" },
            { name:"V.cholerae_07.fa" },
            { name:"V.cholerae_08.fa" },
            { name:"V.cholerae_09.fa" },
            { name:"V.cholerae_010.fa" },
            { name:"V.cholerae_011.fa" },
            { name:"V.cholerae_012.fa" },
        ];

        let i = 0;
        for(i; i < exampleFile.length; i++){
            this.dropzone.emit("addedfile", exampleFile[i]);
            this.dropzone.emit("success", exampleFile[i]);
            this.dropzone.emit("complete", exampleFile[i]);
            this.dropzone.files.push(exampleFile[i]);
        };

        let encodeExampleData = [
            require('./static/Example_data/V.cholerae_01.fa'), 
            require('./static/Example_data/V.cholerae_02.fa'), 
            require('./static/Example_data/V.cholerae_03.fa'), 
            require('./static/Example_data/V.cholerae_04.fa'), 
            require('./static/Example_data/V.cholerae_05.fa'), 
            require('./static/Example_data/V.cholerae_06.fa'), 
            require('./static/Example_data/V.cholerae_07.fa'), 
            require('./static/Example_data/V.cholerae_08.fa'), 
            require('./static/Example_data/V.cholerae_09.fa'), 
            require('./static/Example_data/V.cholerae_10.fa'), 
            require('./static/Example_data/V.cholerae_11.fa'), 
            require('./static/Example_data/V.cholerae_12.fa'), 
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
        let VC07 = new File([decodeExampleData[6]],'V.cholerae_07.fa');
        let VC08 = new File([decodeExampleData[7]],'V.cholerae_08.fa');
        let VC09 = new File([decodeExampleData[8]],'V.cholerae_09.fa');
        let VC10 = new File([decodeExampleData[9]],'V.cholerae_10.fa');
        let VC11 = new File([decodeExampleData[9]],'V.cholerae_11.fa');
        let VC12 = new File([decodeExampleData[9]],'V.cholerae_12.fa');

        let decodeExampleFile = [ VC01, VC02, VC03, VC04, VC05, VC06, VC07, VC08,
            VC09, VC10, VC11, VC12 ];

        let k = 0;
        let upload_status = [];

        for(k; k < decodeExampleFile.length; k++){
            let form = new FormData();
            form.append('file',decodeExampleFile[k]);
            form.append('batch_id',window.batchid);

            fetch('api/profiling/sequence/', {
                method:'POST',
                body:form ,
            }).then(res => upload_status.push(res.status));
        };

        window.databaseName = "Vibrio_cholerae";
        window.tabSwitch = true;

        this.setState(state => ({ switch: true }));

        function trigger_celery(){

            if( upload_status.length == 12 ){
                let scheme = {};
                scheme.occurrence = "95";
                scheme.database = window.databaseName;
                scheme.id = window.batchid;
                fetch('api/profiling/profiling-tree/', {
                    method:'POST',
                    headers: new Headers({'content-type': 'application/json'}),
                    body: JSON.stringify(scheme)
                });
                clearInterval(interval);
            }
        };

        let interval = setInterval(trigger_celery,1500);
    }

    back(){
        window.location.reload();
    }

    
    render() {

        const config = this.componentConfig;
        const djsConfig = this.djsConfig;
        const eventHandlers = {
            init: dz => this.dropzone = dz,}
        const { classes } = this.props;

        if(this.state.profile_result_all == undefined){
            return (
            <div>
                <br />
                <br />
                <div>
                    <div style={{ float:'left' }}>
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
                        <Link to="/profile_view" style={{ textDecoration:'none', color:'#000' }}>
                            Example
                            &nbsp;&nbsp;
                            <CloudUploadIcon />
                        </Link>
                    </Button>
                    &nbsp;&nbsp;&nbsp;&nbsp;
                    <Button variant="contained" className ={classes.cssRoot} 
                     onClick={this.handlePost.bind(this)}>
                        Upload
                        &nbsp;&nbsp;
                        <CloudUploadIcon />
                    </Button>
                    &nbsp;&nbsp;&nbsp;&nbsp;
                    <Link to={this.state.to} style={{ textDecoration:'none' }}>
                        <Button variant="contained" color="primary" onClick={this.submit.bind(this)}>
                            profiling
                            &nbsp;&nbsp;
                            <Icon>send</Icon>
                        </Button>
                    </Link>
                </div>
                <br />
                <br />
                <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <SearchBar
                        onChange = {(value) => this.setState({ querybyID: value })}
                        onRequestSearch={this.query.bind(this)}
                        placeholder = {"Input batch ID here to query result"}
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
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <a download href={this.state.profile_result_all} 
                             style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Download profiles (.tsv)
                                &nbsp;&nbsp;
                                <DownloadIcon />
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
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
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <img src={this.state.svg_file} />
                        </div>
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <font size="4">Download Dendrogram</font>
                        </div>
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            
                            <a download href={this.state.png_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Png 
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.pdf_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Pdf
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.svg_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Svg
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.emf_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                emf
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.newick_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                newick
                                </Button>
                            </a>
                        </div>
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                            <Link to="/" style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default" onClick={this.back.bind(this)}>
                                    <ReplyIcon />
                                    &nbsp;&nbsp;
                                    Back
                                </Button>
                            </Link>
                        </div>
                        <br />
                        <br />
                    </div>
                );
        }
        
    }
}


export default withStyles(styles)(Upload_contigs);