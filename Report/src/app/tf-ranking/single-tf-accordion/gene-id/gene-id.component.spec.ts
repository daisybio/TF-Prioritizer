import { ComponentFixture, TestBed } from '@angular/core/testing';

import { GeneIDComponent } from './gene-id.component';

describe('GeneIDComponent', () => {
  let component: GeneIDComponent;
  let fixture: ComponentFixture<GeneIDComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ GeneIDComponent ]
    })
    .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(GeneIDComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
